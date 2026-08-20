# SA density and seeding backends

minibwa can trade memory for mapping speed along **two independent axes**. This
guide explains both, how to build the indexes for them, and which combination to
pick.

| Axis | Knob | What it changes |
|------|------|-----------------|
| **SA density** | `minibwa index -u INT` | How finely the suffix array is sampled (1/16, 1/8, …). Denser = faster SA lookups, more RAM. |
| **Seeding backend** | `--index-regime …-b2` (needs a `make B2=1` build) | Whether seeding uses minibwa's native BWT or a bwa-mem2/bwa-mem3 **cp_occ** FM-index (an extra SMEM-layer speedup). |

The two are orthogonal: you can run either backend at either density. SA
*density* never changes the alignments — all native densities are byte-identical.
The two *backends* are not guaranteed identical: cp_occ differs from native for
~0.05% of reads (MAPQ-0 multimappers; see §4). Both axes primarily trade lookup
speed for memory.

---

## 1. SA density (`-u`)

The suffix array is sampled at a rate of `1/(1<<u)`. Sparser sampling stores
fewer SA entries (less RAM) but needs more BWT steps per lookup; denser sampling
is the reverse.

- `-u 4` → 1/16 (the historical default; smallest RAM)
- `-u 3` → **1/8 (the current default)** — ~5–6% faster mapping for ~+3 GB RAM
- `-u 2` → 1/4 (denser still; rarely worth it)

The rate is stored in the `.mbw` header, so **existing indexes keep working** —
only newly built indexes change. `-u 4` at map time is not a thing; density is
chosen at **index build time**. (`-u` on `minibwa map` is unrelated — there it
means "suppress unmapped".)

---

## 2. Regimes and the picker

A **regime** is one concrete (density, backend) choice for a given index prefix.
minibwa discovers the regimes available for an index and picks one automatically,
or you can force one.

```
minibwa map --list-regimes  ref reads_1.fq reads_2.fq     # show what's available
minibwa map --index-regime sa8   ref reads_1.fq reads_2.fq # force native 1/8
minibwa map --index-regime sa8-b2 ref reads_1.fq reads_2.fq # force cp_occ 1/8
minibwa map --index-mem 32G ref reads_1.fq reads_2.fq     # auto-pick under a RAM cap
```

Regime names encode the density fraction and backend:

| Name | Backend | Sampling | Notes |
|------|---------|----------|-------|
| `sa16` | native BWT | 1/16 | smallest RAM |
| `sa8`  | native BWT | 1/8  | the default density |
| `sa16-b2` | cp_occ | 1/16 | needs a co-located bwa-mem3 index |
| `sa8-b2`  | cp_occ | 1/8  | needs a co-located bwa-mem3 index |

**Auto-selection** (`--index-regime auto`, the default when you don't pass one)
picks the fastest regime that fits the memory budget. The budget is
`min(host RAM, cgroup limit, --index-mem)`. If nothing fits, it falls back to the
smallest-footprint regime. An explicit `--index-regime NAME` always wins over
auto.

**Mode restriction:** the cp_occ (`-b2`) regimes are **short-read / paired-end
only**. `--meth`, `--hic`, and `--long`/`-x lr` automatically fall back to a
native regime (and forcing `--index-regime sa8-b2 --hic` is rejected).

---

## 3. Building indexes

### Native minibwa index

```
minibwa index ref.fa ref                 # default 1/8
minibwa index -u 4 ref.fa ref            # 1/16 (smaller)
minibwa index -u 3,4 ref.fa ref          # BOTH densities in one pass
```

`-u 3,4` builds the BWT **once** and writes the densest sampling into `ref.mbw`
plus a `ref.sa.u4` sidecar for the other density — so a single build gives you
both `sa8` and `sa16` regimes to choose between at map time, with no duplicated
BWT work.

### cp_occ backend index (optional)

The `-b2` regimes reuse a **bwa-mem2/bwa-mem3 index built at the same prefix**:

```
bwa-mem3 index -p ref ref.fa             # cp_occ FM-index, co-located at prefix "ref"
```

The cp_occ SA sampling rate is fixed by the bwa-mem3 build (there is no `-u` on
`bwa-mem3 index`); it is surfaced in the regime name (`sa8-b2`, `sa16-b2`, …) and
shown by `minibwa map --list-regimes`.

minibwa's `-b2` seeding then loads `ref.bwt.2bit.64`/`.pac`/`.amb`/`.ann`. It
prefers a co-located `ref.l2b` for the reference layer (fast); if none exists it
reconstructs the reference from the bwa-mem3 index itself, so a bwa-mem3-only
index still maps.

The cp_occ backend is **compile-optional** — it is present only in a
`make B2=1` build (linked against a bwa-mem3 `libbwa.a`). A stock `make` build
has no `-b2` regimes at all.

---

## 4. Which should I use?

Measured on hg38, 15M read pairs, 16 threads, warm cache (Apple M-series,
28 cores). Throughput is mapping-phase pairs/s; speedup is versus stock minibwa
(`sa16`, native 1/16) at the same thread count. Ordered fastest first:

| regime | index | sampling | est. RAM | pairs/s | speedup |
|--------|-------|----------|----------|---------|---------|
| `sa1-b2`  | cp_occ | 1/1  | 39 GB | 204,900 | 1.31× |
| `sa2-b2`  | cp_occ | 1/2  | 24 GB | 202,600 | 1.30× |
| `sa4-b2`  | cp_occ | 1/4  | 17 GB | 197,500 | 1.27× |
| `sa1`     | native | 1/1  | 54 GB | 191,400 | 1.23× |
| `sa2`     | native | 1/2  | 30 GB | 189,600 | 1.22× |
| `sa4`     | native | 1/4  | 18 GB | 185,900 | 1.19× |
| `sa8-b2`  | cp_occ | 1/8  | 13 GB | 182,600 | 1.17× |
| `sa16-b2` | cp_occ | 1/16 | 11 GB | 174,900 | 1.12× |
| `sa8`     | native | 1/8  | 12 GB | 174,400 | 1.12× |
| `sa16`    | native | 1/16 |  9 GB | 155,900 | 1.00× |

SA *density* never changes the alignments — it only trades RAM for lookup speed
(all native rows are byte-identical; cp_occ differs from native only in ~0.05%
of reads, MAPQ-0 multimappers, at every density). The speed curve **flattens
sharply past 1/4**: native 1/16 → 1/8 → 1/4 → 1/2 → 1/1 buys
1.00 → 1.12 → 1.19 → 1.22 → 1.23×, so 1/2 and 1/1 cost 2–6× the RAM of 1/8 for
almost nothing (they are included here for completeness, not as recommendations).

Rules of thumb:

- **Most users: take the default (`sa8`, native 1/8).** ~1.12× over 1/16 for
  ~+3 GB, needs only the `minibwa index` you already build, works in a stock
  build.
- **Save memory:** `-u 4` / `--index-regime sa16` (or `--index-mem`) — ~9 GB for
  hg38, ~10% slower than 1/8.
- **A bit more speed at similar RAM:** `sa4` (native 1/4, ~18 GB) is ~1.19×; past
  that the curve is flat, so denser native sampling is not worth the memory.
- **Fastest, if you can spare the bwa-mem3 index + a `make B2=1` build:** the
  cp_occ backend adds ~5–12% over native at the same density (≈5% at 1/8, ≈6% at
  1/4, ≈12% at 1/16). `sa8-b2` is the
  practical pick (~1.17×, 13 GB); `sa4-b2` (~1.27×, 17 GB) is the value sweet
  spot — it beats *native 1/1* (1.23×, 54 GB) at a third of the RAM. Note
  `sa16-b2` (cp_occ 1/16, 11 GB) ≈ `sa8` (native 1/8, 12 GB) in speed — cp_occ
  lets you drop a density notch at lower RAM.

Auto-selection deliberately **does not prefer** a native regime denser than 1/8
(they carry `speed_rank -1`, so the preferred pass skips them); its budget
fallback may still select one when no sparser regime fits. Force a specific
regime with `--index-regime` if you want one of them.

> **Note on cp_occ startup:** if `ref.l2b` is *absent*, the `-b2` backend
> reconstructs the reference from the bwa-mem3 index at load time — an O(genome)
> step that adds several seconds on a full genome. Keep the `.l2b` alongside the
> cp_occ index (the normal case) and startup is ~1s. Startup also grows with
> density (the SA is read at load): ~0.8s at 1/16 up to ~4–5s for a full 1/1 SA —
> negligible for large jobs, but it favours sparser sampling for tiny ones.

---

## Quick recipes

```
# Default, simplest — native 1/8:
minibwa index ref.fa ref
minibwa map -t 16 ref reads_1.fq reads_2.fq > out.sam

# Both densities from one build, then let RAM decide at map time:
minibwa index -u 3,4 ref.fa ref
minibwa map -t 16 --index-mem 10G ref reads_1.fq reads_2.fq > out.sam   # picks sa16

# Fastest short-read PE (needs a `make B2=1` binary):
minibwa index ref.fa ref            # keeps ref.l2b for the fast cp_occ load
bwa-mem3 index -p ref ref.fa
minibwa map -t 16 --index-regime sa8-b2 ref reads_1.fq reads_2.fq > out.sam
```
