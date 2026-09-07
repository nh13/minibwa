<!-- CTA:FEATURES ops-distro,ll-affine-reassoc,ksw2-extension-kernels,single-copy-parser,index-threads-v2,extd2-avx512,inline-appenders,parallel-encode,meth-cleanups,meth-sam-tags,soft-clip-penalty,submem-ablation,sa-default-u3,sa-regime-picker,cp-occ-backend,alt-liftgroup -->
<!-- CTA:MEASURED 2026-09-06 -->
## ⚡ Faster than stock — and ALT-aware

Fully-engaged **nh13/minibwa** vs stock **lh3/minibwa** — same reads, same reference, same compiler:

| workload | arm64 (Graviton4) | x86 (AVX-512) |
|---|---|---|
| short reads (WGS / WES / panel) | **up to 1.23× faster** | up to 1.08× |
| long reads (HiFi / ONT) | up to 1.16× | **up to 1.26× faster** |
| Hi-C | up to 1.09× | up to 1.02× |
| ALT-shadowed reads placed confidently | **6.1× more than stock** | **6.1× more than stock** |

<sub>Whole-distribution defaults vs stock lh3/minibwa · hs38DH · <code>-t 16</code>, 2 reps · lh3 has no ALT handling · output byte-identical to upstream except where a feature is engaged (≥99.95% concordance) · measured 2026-09-06.</sub>
