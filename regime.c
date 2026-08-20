#include "regime.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <glob.h>
#if defined(__APPLE__)
#include <sys/sysctl.h>
#endif

/* Fixed safety margin added on top of the raw on-disk footprint when
 * estimating a regime's peak RAM: covers query buffers, thread-local
 * scratch, and the fact that l2b/mbw file sizes are a conservative proxy
 * for their in-memory representation (over-estimating is the safe
 * direction here). */
#define MB_REGIME_MARGIN 2684354560ULL /* 2.5 GiB */

int mb_regime_pick(const mb_regime_t *r, int n, uint64_t budget, uint32_t mode, const char *forced){
	int best = -1;
	for (int i = 0; i < n; ++i) {
		if (!(r[i].mode_mask & mode)) continue;
		if (forced) { if (strcmp(r[i].name, forced)==0) return i; else continue; }
		if (r[i].est_ram > budget) continue;
		if (best < 0 || r[i].speed_rank > r[best].speed_rank) best = i;
	}
	return best;
}

static uint64_t host_avail_bytes(void){
#if defined(__APPLE__)
	int mib[2] = { CTL_HW, HW_MEMSIZE }; uint64_t v = 0; size_t len = sizeof v;
	if (sysctl(mib, 2, &v, &len, NULL, 0) == 0) return v;
	return 0;
#else
	FILE *fp = fopen("/proc/meminfo","r"); char k[64]; unsigned long kb;
	if (!fp) return 0;
	while (fscanf(fp,"%63s %lu kB\n",k,&kb)==2)
		if (strcmp(k,"MemAvailable:")==0){ fclose(fp); return (uint64_t)kb*1024ULL; }
	fclose(fp); return 0;
#endif
}
static uint64_t cgroup_limit_bytes(void){
	/* cgroup v2 then v1; "max"/absurd values => no limit (0) */
	unsigned long long v; FILE *fp;
	if ((fp=fopen("/sys/fs/cgroup/memory.max","r"))){
		char b[64]; if (fgets(b,sizeof b,fp)){ fclose(fp);
			if (!strncmp(b,"max",3)) return 0;
			v=strtoull(b,0,10); return (v && v < (1ULL<<62))? v : 0; }
		fclose(fp);
	}
	if ((fp=fopen("/sys/fs/cgroup/memory/memory.limit_in_bytes","r"))){
		if (fscanf(fp,"%llu",&v)==1){ fclose(fp); return (v && v < (1ULL<<62))? v : 0; }
		fclose(fp);
	}
	return 0;
}
uint64_t mb_mem_budget(uint64_t user_cap){
	uint64_t host = host_avail_bytes(), cg = cgroup_limit_bytes();
	uint64_t b = host;
	if (cg && cg < b) b = cg;
	if (user_cap && user_cap < b) b = user_cap;
	return b;
}

/* Read a big/little-endian-agnostic (host-native) uint32 at byte offset
 * `offset` in `path`. Returns 0 on success, -1 if the file can't be opened,
 * seeked, or doesn't have enough bytes there. */
static int read_u32_at(const char *path, long offset, uint32_t *out){
	FILE *fp = fopen(path, "rb");
	if (!fp) return -1;
	if (fseek(fp, offset, SEEK_SET) != 0) { fclose(fp); return -1; }
	size_t n = fread(out, 4, 1, fp);
	fclose(fp);
	return n == 1 ? 0 : -1;
}

/* speed_rank for a BWT-backend regime: denser sampling (smaller sa_bit) is
 * faster to query, so it ranks higher. sa_bit < 3 (denser than 1/8) is
 * off the supported frontier for M2 -- such a regime may still be present
 * (e.g. a leftover sidecar) but must never be auto-picked, hence rank -1. */
static int bwt_speed_rank(int sa_bit){
	return sa_bit < 3 ? -1 : 2 * (6 - sa_bit);
}

static void fill_bwt_regime(mb_regime_t *rg, int sa_bit, uint64_t est_ram, const char *sa_path){
	memset(rg, 0, sizeof *rg);
	rg->backend = MB_BACKEND_BWT;
	rg->sa_bit = sa_bit;
	rg->est_ram = est_ram;
	rg->speed_rank = bwt_speed_rank(sa_bit);
	rg->mode_mask = MB_MODE_SRPE | MB_MODE_METH | MB_MODE_HIC | MB_MODE_LR;
	snprintf(rg->name, sizeof rg->name, "sa%d", 1 << sa_bit);
	if (sa_path) strncpy(rg->sa_path, sa_path, sizeof rg->sa_path - 1);
	else rg->sa_path[0] = '\0';
}

/* Discover on-disk SA-density regimes for `prefix`:
 *  - the density bundled in <prefix>.mbw (sa_path == ""), read from the
 *    4-byte sa_bit field at offset 4 of the .mbw header; and
 *  - one regime per <prefix>.sa.u* sidecar, read the same way from offset 4
 *    of the sidecar (its header is magic[4], sa_bit u32 @4, n_sa u64 @8, ...).
 * A sidecar whose sa_bit matches the bundled density is skipped -- it is a
 * redundant duplicate of the bundled regime, so listing it again would just
 * clutter `mb_regime_list_print` and `auto` picking without adding a
 * reachable option.
 * b2_available is accepted for forward compatibility with M3's cp_occ
 * regimes; nothing is discovered for it in M2. */
int mb_regime_discover(const char *prefix, int b2_available, mb_regime_t *out, int max){
	char fn_l2b[1152], fn_mbw[1152];
	struct stat st_l2b, st_mbw;
	uint32_t bundled_sa_bit;
	uint64_t l2b_size, mbw_size;
	int n = 0;

	if (!prefix || !out || max <= 0) return 0;

	snprintf(fn_l2b, sizeof fn_l2b, "%s.l2b", prefix);
	snprintf(fn_mbw, sizeof fn_mbw, "%s.mbw", prefix);

	if (stat(fn_l2b, &st_l2b) != 0 || stat(fn_mbw, &st_mbw) != 0) return 0;
	if (read_u32_at(fn_mbw, 4, &bundled_sa_bit) != 0) return 0;

	l2b_size = (uint64_t)st_l2b.st_size;
	mbw_size = (uint64_t)st_mbw.st_size;

	/* The bundled regime: its SA lives inside .mbw itself. */
	if (n < max) {
		fill_bwt_regime(&out[n], (int)bundled_sa_bit, l2b_size + mbw_size + MB_REGIME_MARGIN, NULL);
		n++;
	}

	/* Sidecar regimes: <prefix>.sa.u* */
	{
		char pattern[1152];
		glob_t gl;
		snprintf(pattern, sizeof pattern, "%s.sa.u*", prefix);
		memset(&gl, 0, sizeof gl);
		if (glob(pattern, 0, NULL, &gl) == 0) {
			size_t i;
			for (i = 0; i < gl.gl_pathc && n < max; ++i) {
				const char *side = gl.gl_pathv[i];
				struct stat st_side;
				uint32_t side_sa_bit;
				if (stat(side, &st_side) != 0) continue;
				if (read_u32_at(side, 4, &side_sa_bit) != 0) continue;
				if (side_sa_bit == bundled_sa_bit) continue; /* de-dup: redundant with the bundled regime */
				fill_bwt_regime(&out[n], (int)side_sa_bit,
					l2b_size + mbw_size + (uint64_t)st_side.st_size + MB_REGIME_MARGIN, side);
				n++;
			}
		}
		globfree(&gl);
	}

	if (b2_available) {
		/* M3 adds cp_occ regimes here (a second SA-lookup backend that
		 * trades RAM for speed via a checkpointed OCC structure). Not
		 * implemented in M2: b2_available is always 0 at call sites. */
	}

	return n;
}

void mb_regime_list_print(FILE *fp, const mb_regime_t *r, int n){
	int i;
	fprintf(fp, "%-10s %-8s %-10s %12s %6s  %s\n",
		"name", "backend", "sampling", "est-RAM(GB)", "rank", "modes");
	for (i = 0; i < n; ++i) {
		char modes[32];
		size_t len;
		modes[0] = '\0';
		if (r[i].mode_mask & MB_MODE_SRPE) strcat(modes, "srpe,");
		if (r[i].mode_mask & MB_MODE_METH) strcat(modes, "meth,");
		if (r[i].mode_mask & MB_MODE_HIC)  strcat(modes, "hic,");
		if (r[i].mode_mask & MB_MODE_LR)   strcat(modes, "lr,");
		len = strlen(modes);
		if (len > 0) modes[len-1] = '\0'; /* trim the trailing comma */
		fprintf(fp, "%-10s %-8s 1/%-8d %12.2f %6d  %s\n",
			r[i].name,
			r[i].backend == MB_BACKEND_BWT ? "bwt" : "cp_occ",
			1 << r[i].sa_bit,
			(double)r[i].est_ram / (double)(1ULL<<30),
			r[i].speed_rank,
			modes);
	}
}
