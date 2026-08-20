/* b2idx — minibwa seeding on the bwa-mem2/bwa-mem3 FM-index. See b2idx.h.
 *
 * This is a faithful port of minibwa's FM-index seeding leaves (mb_bwt_extend,
 * mb_bwt_set_intv, mb_bwt_count_kmer, mb_bwt_smem_batch, mb_bwt_sa_batch) with
 * the single occ primitive swapped from minibwa's classic-BWT rank to
 * bwa-mem2's checkpointed cp_occ. The algorithm, the SMEM state machine, the
 * 10-mer cache, and the cross-read block prefetch are structurally identical to
 * minibwa's; only the leaf occ/SA reads hit the bwa-mem2 index. Lifted from the
 * validated seed-bench spike (mb_on_b3b: 100% batch-neutral vs the scalar path,
 * 99.25% seed-identical to stock minibwa, the 0.75% being benign N-substitution
 * index-content differences). */

#include "b2idx.h"          /* mb_sai_t, mb_sai_v, mb_smem_entry_t + the C API */
#include "kalloc.h"         /* Kgrow / kmalloc / kfree (minibwa's pool) */
#include <vector>
#include <cstring>
#include <cstdint>
#include <cstdio>

#include "fmi_seed_api.h"   /* FmiSeed, SMEM, CP_OCC, one_hot_mask_array,
                               CP_SHIFT/CP_MASK, _mm_countbits_64,
                               fmi_seed_open/close/cp_occ/count/sentinel/sa_prefetch */

/* fmi_seed_api.h is a lean facade and does not pull in NEON intrinsics (unlike
 * FMI_search.h, which drags simd_compat.h in transitively). b2_occ_sp_sz below
 * uses raw NEON types/intrinsics on arm64, so include them directly here. */
#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
#include <arm_neon.h>
#endif

/* ------------------------------- handle ----------------------------------- */

struct mb_b2_t {
    FmiSeed       *fmi;
    const CP_OCC  *cp_occ;
    const int64_t *count;     /* C() array, sentinel already folded in */
    int64_t        sentinel;  /* S^{-1}(0) */
    std::vector<mb_sai_t> pre; /* 10-mer cache: pre[kmer] -> interval */
    int            pre_len;
};

/* ------------------------- cp_occ FMD primitives -------------------------- */

/* occ_sp[b] (occ at sp) and sz[b] (occ at ep minus occ at sp) for all four bases,
 * from the two checkpoint blocks. Bit-identical to bwa-mem2's backwardExt: on
 * arm64 it uses the same NEON in-lane path (all four bases' popcount in-vector,
 * no GPR<->SIMD round-trips — the round-trips are why the arm64 scalar popcount
 * loses); elsewhere the scalar hardware-popcount loop, which is already optimal
 * on x86. This is the occ hot leaf, called ~10^9x per WGS run. */
static inline void b2_occ_sp_sz(const CP_OCC *cp_occ, int64_t sp, int64_t ep,
                                int64_t occ_sp[4], int64_t sz[4])
{
    const CP_OCC &blk_sp = cp_occ[sp >> CP_SHIFT];
    const CP_OCC &blk_ep = cp_occ[ep >> CP_SHIFT];
    const uint64_t mask_sp = one_hot_mask_array[sp & CP_MASK];
    const uint64_t mask_ep = one_hot_mask_array[ep & CP_MASK];
#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
    const uint64x2_t msp = vdupq_n_u64(mask_sp);
    const uint64x2_t mep = vdupq_n_u64(mask_ep);
    #define B2_PC64(v) vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(vcntq_u8(vreinterpretq_u8_u64(v)))))
    const uint64x2_t psp01 = B2_PC64(vandq_u64(vld1q_u64(&blk_sp.one_hot_bwt_str[0]), msp));
    const uint64x2_t psp23 = B2_PC64(vandq_u64(vld1q_u64(&blk_sp.one_hot_bwt_str[2]), msp));
    const uint64x2_t pep01 = B2_PC64(vandq_u64(vld1q_u64(&blk_ep.one_hot_bwt_str[0]), mep));
    const uint64x2_t pep23 = B2_PC64(vandq_u64(vld1q_u64(&blk_ep.one_hot_bwt_str[2]), mep));
    #undef B2_PC64
    const uint64x2_t occ_sp01 = vaddq_u64(vld1q_u64((const uint64_t*)&blk_sp.cp_count[0]), psp01);
    const uint64x2_t occ_sp23 = vaddq_u64(vld1q_u64((const uint64_t*)&blk_sp.cp_count[2]), psp23);
    const uint64x2_t occ_ep01 = vaddq_u64(vld1q_u64((const uint64_t*)&blk_ep.cp_count[0]), pep01);
    const uint64x2_t occ_ep23 = vaddq_u64(vld1q_u64((const uint64_t*)&blk_ep.cp_count[2]), pep23);
    vst1q_u64((uint64_t*)&occ_sp[0], occ_sp01);
    vst1q_u64((uint64_t*)&occ_sp[2], occ_sp23);
    vst1q_u64((uint64_t*)&sz[0], vsubq_u64(occ_ep01, occ_sp01));
    vst1q_u64((uint64_t*)&sz[2], vsubq_u64(occ_ep23, occ_sp23));
#else
    for (int b = 0; b < 4; b++) {
        int64_t o_s = blk_sp.cp_count[b] + _mm_countbits_64(blk_sp.one_hot_bwt_str[b] & mask_sp);
        int64_t o_e = blk_ep.cp_count[b] + _mm_countbits_64(blk_ep.one_hot_bwt_str[b] & mask_ep);
        occ_sp[b] = o_s; sz[b] = o_e - o_s;
    }
#endif
}

/* mb_bwt_set_intv over bwa-mem2 count[] (== minibwa L2[c]+1). */
static inline void b2_set_intv(const mb_b2_t *b2, int c, mb_sai_t *ik)
{
    ik->x[0] = (uint64_t)b2->count[c];
    ik->x[1] = (uint64_t)b2->count[3 - c];
    ik->size = (uint64_t)(b2->count[c + 1] - b2->count[c]);
    ik->info = 0;
}

/* mb_bwt_extend re-expressed over cp_occ; fills all four next intervals.
 * Ranked coord = ik->x[!is_back] -> count[b]+occ; cumulated coord = ik->x[is_back]
 * accumulates from base 3 down with the sentinel fix folded into base 3. */
static inline void b2_extend(const mb_b2_t *b2, const mb_sai_t *ik, mb_sai_t ok[4], int is_back)
{
    const int fwd = !is_back, rev = is_back;
    const int64_t base = (int64_t)ik->x[fwd];
    const int64_t s0   = (int64_t)ik->size;
    int64_t occ_sp[4], sz[4];
    b2_occ_sp_sz(b2->cp_occ, base, base + s0, occ_sp, sz);
    for (int b = 0; b < 4; b++) ok[b].x[fwd] = (uint64_t)(b2->count[b] + occ_sp[b]);
    const int64_t soff = (base <= b2->sentinel && base + s0 > b2->sentinel) ? 1 : 0;
    ok[3].x[rev] = ik->x[rev] + (uint64_t)soff;
    ok[2].x[rev] = ok[3].x[rev] + (uint64_t)sz[3];
    ok[1].x[rev] = ok[2].x[rev] + (uint64_t)sz[2];
    ok[0].x[rev] = ok[1].x[rev] + (uint64_t)sz[1];
    for (int b = 0; b < 4; b++) ok[b].size = (uint64_t)sz[b];
}

/* Prefetch the checkpoint block covering SA position k (one CP_OCC per line). */
static inline void b2_prefetch(const mb_b2_t *b2, uint64_t k)
{
    __builtin_prefetch(&b2->cp_occ[k >> CP_SHIFT]);
}

/* ---------------------------- 10-mer cache -------------------------------- */
/* Port of mb_bwt_count_kmer (bwt.c). The kmer bit convention here is paired
 * with the batch lookup below, exactly as in minibwa. */
struct b2_kstack_t { mb_sai_t p; int32_t d; uint8_t c; };

static void b2_build_cache(mb_b2_t *b2, int depth)
{
    b2->pre_len = depth;
    mb_sai_t zero; std::memset(&zero, 0, sizeof(zero));
    b2->pre.assign((size_t)1 << (depth * 2), zero);
    mb_sai_t *s = b2->pre.data();
    b2_kstack_t stack[64];
    uint8_t str[16];
    int s_top = 0;
    for (int a = 0; a < 4; ++a) {
        b2_kstack_t *p = &stack[s_top++];
        b2_set_intv(b2, a, &p->p);
        p->d = 1; p->c = (uint8_t)a;
    }
    while (s_top > 0) {
        b2_kstack_t top = stack[--s_top];
        mb_sai_t ok[4];
        if (top.d > 0) str[depth - top.d] = top.c;
        b2_extend(b2, &top.p, ok, 1);
        for (int a = 0; a < 4; ++a) {
            str[depth - top.d - 1] = (uint8_t)a;
            if (top.d != depth - 1) {
                b2_kstack_t *p = &stack[s_top++];
                p->p = ok[a]; p->d = top.d + 1; p->c = (uint8_t)a;
            } else {
                uint64_t x = 0;
                for (int i = 0; i < depth; ++i) x |= (uint64_t)str[i] << (i * 2);
                s[x] = ok[a];
            }
        }
    }
}

/* --------------------------- SMEM state machine --------------------------- */
/* Exact structural port of mb_bwt_smem_batch (bwt.c:356), operating on the
 * caller's mb_smem_entry_t array (honoring per-entry min_len/min_occ/st/en) and
 * appending to each entry's v via the caller's kalloc pool. */

static inline void b2_one_step_back(const mb_b2_t *b2, mb_smem_entry_t *s)
{
    mb_sai_t ok[4];
    int32_t c = s->q[s->i];
    b2_extend(b2, &s->p, ok, 1);
    if (ok[c].size < (uint64_t)s->min_occ) { s->x = s->i + 1; s->stage = 1; }
    else {
        s->p = ok[c]; s->i--;
        b2_prefetch(b2, s->p.x[0]);
        b2_prefetch(b2, s->p.x[0] + s->p.size);
    }
}

extern "C" void mb_b2_smem_batch(void *km, void *b2_, int32_t n, mb_smem_entry_t *a)
{
    const mb_b2_t *b2 = (const mb_b2_t *)b2_;
    const mb_sai_t *pre = b2->pre.data();
    const int pre_len = b2->pre_len;

    /* init: stage 1, x=st, preallocate each v to >=64 (mirrors minibwa).
     * Reuse a thread-local ring queue so the hot path allocates nothing. */
    static thread_local std::vector<int32_t> q;
    q.clear(); q.reserve(n);
    for (int32_t i = 0; i < n; ++i) {
        mb_smem_entry_t *s = &a[i];
        s->stage = 1;
        s->x = s->st;
        if (s->v->m < 64) { s->v->m = 64; s->v->a = (mb_sai_t*)krealloc(km, s->v->a, s->v->m * sizeof(mb_sai_t)); }
        q.push_back(i);
    }

    size_t head = 0;
    while (head < q.size()) {
        int32_t idx = q[head++];
        mb_smem_entry_t *s = &a[idx];
        if (s->stage == 1) {
            int32_t i, xn;
            if (s->en - s->x < s->min_len) continue;   /* drop: skip re-queue */
            for (i = s->x, xn = -1; i < s->x + s->min_len; ++i)
                if (s->q[i] > 3) xn = i;
            if (xn >= 0) { s->x = xn + 1; }
            else {
                s->i = s->x + s->min_len - 1;
                if (pre_len && s->min_len >= pre_len) {
                    for (i = 0, s->kmer = 0; i < pre_len; ++i, s->i--)
                        s->kmer = s->kmer << 2 | s->q[s->i];
                    __builtin_prefetch(&pre[s->kmer]);
                    s->stage = 2;
                } else {
                    b2_set_intv(b2, s->q[s->i--], &s->p);
                    s->stage = 3;
                }
            }
        } else if (s->stage == 2 || s->stage == 5) {
            s->p = pre[s->kmer];
            if (s->p.size < (uint64_t)s->min_occ) {
                s->i += pre_len;
                b2_set_intv(b2, s->q[s->i--], &s->p);
            }
            b2_prefetch(b2, s->p.x[0]);
            b2_prefetch(b2, s->p.x[0] + s->p.size);
            s->stage++;
        } else if (s->stage == 3) {
            if (s->i < s->x) {
                b2_prefetch(b2, s->p.x[1]);
                b2_prefetch(b2, s->p.x[1] + s->p.size);
                s->i = s->x + s->min_len;
                s->stage = 4;
            } else b2_one_step_back(b2, s);
        } else if (s->stage == 4) {
            if (s->i == s->en) {
                s->p.info = (uint64_t)s->x << 32 | (uint32_t)s->i;
                Kgrow(km, mb_sai_t, s->v->a, s->v->n, s->v->m);
                s->v->a[s->v->n++] = s->p;
                continue;
            } else {
                int32_t i, c = 3 - (int32_t)s->q[s->i];
                mb_sai_t ok[4];
                if (c >= 0) b2_extend(b2, &s->p, ok, 0);
                if (c >= 0 && ok[c].size >= (uint64_t)s->min_occ) {
                    s->p = ok[c]; s->i++;
                    b2_prefetch(b2, s->p.x[1]);
                    b2_prefetch(b2, s->p.x[1] + s->p.size);
                } else {
                    s->p.info = (uint64_t)s->x << 32 | (uint32_t)s->i;
                    Kgrow(km, mb_sai_t, s->v->a, s->v->n, s->v->m);
                    s->v->a[s->v->n++] = s->p;
                    if (c < 0) { s->x = s->i + 1; s->stage = 1; }
                    else if (pre_len && s->i - s->x - 1 >= pre_len) {
                        for (i = 0, s->kmer = 0; i < pre_len; ++i, s->i--)
                            s->kmer = s->kmer << 2 | s->q[s->i];
                        __builtin_prefetch(&pre[s->kmer]);
                        s->stage = 5;
                    } else {
                        b2_set_intv(b2, s->q[s->i--], &s->p);
                        s->stage = 6;
                    }
                }
            }
        } else if (s->stage == 6) {
            if (s->i < s->x + 1) { s->x = s->i + 1; s->stage = 1; }
            else b2_one_step_back(b2, s);
        }
        q.push_back(idx);
        if (head > (1u << 16) && head * 2 > q.size()) {
            q.erase(q.begin(), q.begin() + head);
            head = 0;
        }
    }
}

/* ------------------------------ SA lookup --------------------------------- */
/* Resolve SA indices to 2N positions via bwa-mem2's prefetch-pooled path. Each
 * a[i] becomes a unit SMEM (k=a[i], s=1); step==1 so no decimation, and the
 * coordArray comes back in input order. */
extern "C" void mb_b2_sa_batch(void *km, void *b2_, int64_t n, uint64_t *a)
{
    (void)km;
    mb_b2_t *b2 = (mb_b2_t *)b2_;
    if (n <= 0) return;
    /* Thread-local reusable staging buffers: this is called once per ~20-coord
     * process_batch, so per-call heap churn would dominate. */
    static thread_local std::vector<SMEM>    sm;
    static thread_local std::vector<int64_t> coord;
    if ((int64_t)sm.size()    < n) sm.resize(n);
    if ((int64_t)coord.size() < n) coord.resize(n);
    for (int64_t i = 0; i < n; i++) {
        std::memset(&sm[i], 0, sizeof(SMEM));
        sm[i].k = (int64_t)a[i];
        sm[i].s = 1;
    }
    int64_t dummy = 0, id = 0;
    fmi_seed_sa_prefetch(b2->fmi, sm.data(), coord.data(), &dummy, n,
                         /*max_occ=*/1, /*tid=*/0, &id);
    for (int64_t i = 0; i < n; i++) a[i] = (uint64_t)coord[i];
}

/* ------------------------------ load / free ------------------------------- */

extern "C" void *mb_b2_load(const char *prefix)
{
    mb_b2_t *b2 = new mb_b2_t();
    b2->fmi = fmi_seed_open(prefix);
    b2->cp_occ   = fmi_seed_cp_occ(b2->fmi);
    b2->count    = fmi_seed_count(b2->fmi);
    b2->sentinel = fmi_seed_sentinel(b2->fmi);
    b2_build_cache(b2, 10);   /* match minibwa's mb_bwt_cache(bwt, 10) */
    return b2;
}

extern "C" void mb_b2_destroy(void *b2_)
{
    mb_b2_t *b2 = (mb_b2_t *)b2_;
    if (!b2) return;
    fmi_seed_close(b2->fmi);   /* frees cp_occ / sa_ms_byte / sa_ls_word */
    delete b2;
}
