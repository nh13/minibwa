
<!-- DISTRO:BEGIN -->
> **This is a downstream build of [lh3/minibwa](https://github.com/lh3/minibwa).**
> It carries changes upstream declined, plus a few not yet offered. Upstream remains the source
> of truth for everything else.
>
> **Install from a release tag, not from a branch.** Tags such as `v0.6-nh13.1` are immutable;
> the `dist` branch is rebuilt and force-pushed on every upstream change and will rewind under
> you. `minibwa version` reports the downstream version, and it appears in the `@PG VN:` tag of
> every SAM file this build writes, so output is always traceable to the build that produced it.
>
> Bug reports for anything in the table below belong here, not upstream.
> Changes investigated and deliberately not shipped are in [`GRAVEYARD.md`](GRAVEYARD.md).

| feature | upstream | status | output | summary |
|---|---|---|---|---|
| `ops-distro` | — | unsubmitted | identical | distribution manifest, assembly engine and workflows |
| `ll-affine-reassoc` | [lh3/minibwa#37](https://github.com/lh3/minibwa/pull/37) | rejected | identical | reassociate the affine-gap recurrence in ksw2_ll for arm64 |
| `ksw2-extension-kernels` | [lh3/minibwa#64](https://github.com/lh3/minibwa/pull/64) | open | identical | shuffle-LUT prepass, one-vext rail shift, 16-wide exact-max and NEON BIT direction bytes in ksw2_extd2/extz2; +11% HiFi and +10% ONT on arm64, +7% on x86 |
| `single-copy-parser` | — | unsubmitted | identical | read FASTQ records with a single copy instead of two; 0.67-1.24 pp whole-stack ladder increment |
| `index-threads-v2` | — | unsubmitted | identical | parallelize the SA-to-BWT pipeline (OpenMP, gated on the existing LIBSAIS_OPENMP probe) |
| `extd2-avx512` | [lh3/minibwa#20](https://github.com/lh3/minibwa/pull/20) | open | identical | AVX2/AVX-512 ksw_extd2 with runtime dispatch; pays on HiFi/ONT |
| `inline-appenders` | [lh3/minibwa#32](https://github.com/lh3/minibwa/pull/32) | rejected | identical | format SAM/PAF records with inline appenders; with parallel-encode worth +3.7 pp at the default -K on wgs-1M but ~0 above ~10M pairs (input-size dependent, see notes) |
| `parallel-encode` | [lh3/minibwa#33](https://github.com/lh3/minibwa/pull/33) | open | identical | format SAM/PAF in the mapping step instead of the output thread; relocates work rather than removing it, so it pays only on inputs small enough to be starvation-bound |
| `meth-cleanups` | [lh3/minibwa#19](https://github.com/lh3/minibwa/pull/19) | rejected | identical | a --meth CI test; documents the mb_align1_inv strand flip and b_ts (comments only) |
| `meth-sam-tags` | [lh3/minibwa#15](https://github.com/lh3/minibwa/pull/15) | open | conditional | emit Bismark-compatible XR/XG/XM tags |
| `soft-clip-penalty` | [lh3/minibwa#13](https://github.com/lh3/minibwa/pull/13) | rejected | conditional | 5'/3' soft-clip penalty (-L) |
| `submem-ablation` | [lh3/minibwa#12](https://github.com/lh3/minibwa/pull/12) | rejected | conditional | expose --max-sub-occ and --min-sub-occ ablation flags |
| `alt-liftgroup` | — | unsubmitted | conditional | ALT-aware mapping via post-extension liftover groups |
<!-- DISTRO:END -->

[![GitHub Downloads](https://img.shields.io/github/downloads/lh3/minibwa/total.svg?style=social&logo=github&label=Download)](https://github.com/lh3/minibwa/releases)
[![Bioconda](https://img.shields.io/conda/dn/bioconda/minibwa.svg?style=flag&label=bioconda)](https://bioconda.github.io/recipes/minibwa/README.html)
[![Homebrew](https://img.shields.io/homebrew/v/minibwa)](https://formulae.brew.sh/formula/minibwa)
[![Build Status](https://github.com/lh3/minibwa/actions/workflows/build.yml/badge.svg)](https://github.com/lh3/minibwa/actions)
[![preprint](https://img.shields.io/badge/arXiv-2606.15357-blue)](https://arxiv.org/abs/2606.15357)

## Getting Started
```sh
git clone https://github.com/lh3/minibwa
cd minibwa && make

# with test data
./minibwa index test/chrM-human.fa.gz chrM-human              # index the genome
./minibwa map chrM-human test/chrM-read_?.fa.gz > aln.sam     # align and output in SAM

# other examples without test data
minibwa map -ft16 ref.index long-read.fq > aln.paf            # align long reads
minibwa map --hic ref.index reads.interleaved.fq > aln.sam    # align Hi-C short reads

# align *directional* bisulfite sequencing (BS-seq) reads
minibwa index --meth -t8 ref.fa                               # generate BS-seq index
minibwa map --meth ref.fa read1.fq read2.fq > aln.sam         # map BS-seq reads
```

## Introduction

Minibwa aligns short reads against a reference genome. It is the successor of
[bwa-mem][bwa] with a different algorithm. Minibwa is over three times as fast as the
original bwa-mem and twice as fast as [bwa-mem2][bwa-mem2] at comparable accuracy. While
minibwa works with accurate long reads, [minimap2][mm2] is more robust under high
error rate.

Minibwa is a hybrid of bwa-mem and minimap2: it indexes the genome with
Burrow-Wheeler Transform (BWT), finds variable-length seeds like bwa-mem, and
performs chaining and SIMD-based nucleotide alignment with the minimap2
algorithm. Minibwa speeds up bwa-mem2 further with additional prefetch for
seeding, new heuristics to skip unnecessary mate rescue and reduced effort in
highly repetitive regions where reads would often be wrongly mapped due to
structural changes anyway.

## Users' Guide

### Intended use cases

Minibwa is designed for mapping short reads and accurate long reads. It does
not support spliced alignment and has not been tuned for aligning long contigs.
For now, minibwa does not properly work with alternate contigs in the reference
genome. Please use a version of the reference without such contigs.

### Installation

Minibwa requires either SSE4.2 on x86 CPUs or NEON on ARM. It depends on
[zlib][zlib] installed on your system and also includes slightly modified
source code of [mimalloc][mimalloc] and [libsais][libsais] which optionally
uses OpenMP for multi-threading. You can build minibwa with
```sh
make             # automatically detect OpenMP and arm64 vs. x86_64
make omp=0       # disable multi-threading in libsais (no effect on mapping)
make gpl=0       # disable GPL'd code for low-memory BWT building (no effect on mapping)
make mimalloc=0  # disable mimalloc and use the system malloc+kalloc instead
```
This produces a single binary `minibwa` which you can copy to your `PATH`.

### Usage

Like bwa-mem, minibwa requires to index the genome before read alignment.

#### Indexing

You can index the reference genome with
```sh
minibwa index -t8 ref.fa     # index with 8 threads, using 18N RAM (N is the genome size)
minibwa index ref.fa prefix  # use a different index prefix instead of ref.fa
minibwa index -l ref.fa      # use less memory at the cost of performance
minibwa index --meth ref.fa  # generate BS-seq index
```
Minibwa generates two files: `ref.fa.l2b` for 2-bit encoded reference genome
sequences and `ref.fa.mbw` for BWT and sampled suffix array. In the `--meth`
mode, minibwa additionally generates `ref.fa.meth.mbw` for the BWT of the
3-base genome.

#### Mapping

By default, minibwa dynamically changes multiple internal parameters based on
individual read lengths. It works for both short and accurate long reads.
```sh
minibwa map -t8 ref.fa read1.fq read2.fq    # map paired-end reads and output SAM
minibwa map -ft8 ref.fa read.fa.gz          # map single-end or long reads; output PAF
minibwa map --hic ref.fa hic1.fq hic2.fq    # map Hi-C short reads
minibwa map --meth ref.fa read1.fq read2.fq # map BS-seq reads; requiring "index --meth"
```
Note in the default adaptive mode, `-g`/`-w`/`-W`/`-N`/`-m`/`-s` only changes
the short-read setting; the long-read setting is fixed. This mode is disabled
with `--adap=no` or when `-x sr` or `-x lr` is specified.

#### Mapping with legacy bwa-mem CLI

Minibwa also provides legacy bwa-mem command-line interface (CLI) via the `mem` subcommand.
However, due to algorithm and parameter differences, many bwa-mem options are ignored.
The output minibwa alignment is also not identical to bwa-mem.

## Developers' Guide

Minibwa provides basic APIs for loading index and aligning reads.
[api-test/ex-one.c](api-test/ex-one.c) shows an example to align each read
independently; [api-test/ex-batch.c](api-test/ex-batch.c) aligns multiple reads
in batch, which is faster and also supports paired-end mapping.
[dev.md](dev.md) explains how minibwa differs from BWA-MEM and minimap2.

## License

Minibwa is distributed under the MIT license. It also incorporates source code
from the following projects:

 * libsais: Apache 2 License. Copyright (c) 2021-2025 Ilya Grebnov
 * mimalloc: MIT License. Copyright (c) 2018-2026 Microsoft Corporation, Daan Leijen

The master branch is optionally built on the following projects:

 * QSufSort: HPND License. Copyright (c) 1999 N. Jesper Larsson
 * bwtgen: GPL 2 License. Copyright (c) 2004 Wong Chi Kwong

Notably, the master branch includes GPL'd [bwtgen.c](bwtgen.c) for low-memory
BWT construction. If you compile this file, which is the default, the resulting
binary will be GPL'd. You can disable the low-memory algorithm with `make
gpl=0` to generate non-GPL binary. The [Apache2 branch][apache2] does not
include GPL'd source code.

## Limitations

* Minibwa does not work with noisy long reads or spliced RNA-seq reads.
* Minibwa does not support undirectional bisulfite sequencing data.
* Minibwa does not recognize alternate haplotypes.

[apache2]: https://github.com/lh3/minibwa/tree/Apache2
[zlib]: https://zlib.net/
[mimalloc]: https://github.com/microsoft/mimalloc
[libsais]: https://github.com/IlyaGrebnov/libsais
[bwa]: https://github.com/lh3/bwa
[mm2]: https://github.com/lh3/minimap2
[bwa-mem2]: https://github.com/bwa-mem2/bwa-mem2
