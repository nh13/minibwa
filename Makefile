CC=			gcc
CFLAGS=		-std=c99 -g -Wall -O3
CPPFLAGS=
LDFLAGS=
INCLUDES=
LOBJS=		kommon.o kalloc.o bwt.o l2bit.o options.o seed.o map-algo.o lchain.o align.o pe.o cs.o format.o \
			ksw2_extz2_sse.o ksw2_extd2_sse.o ksw2_ll_sse.o regime.o
AOBJS=		kthread.o libsais.o libsais64.o index.o bseq.o map-main.o fastmap.o
MALLOC_O=	mimalloc.o
PROG=		minibwa
LIBS=		-lpthread -lz -lm
ARCH=		$(shell uname -m)
omp=		$(shell printf '\043include <omp.h>\nint main(){return 0;}' | $(CC) -x c -fopenmp -o /dev/null - 2>/dev/null && echo "1" || echo "0")

ifneq ($(asan),)
	override CFLAGS+=-fsanitize=address
	override LDFLAGS+=-fsanitize=address
	override LIBS+=-ldl
endif

ifeq ($(omp),1)
	override CPPFLAGS+=-DLIBSAIS_OPENMP
	override CFLAGS+=-fopenmp
	override LIBS+=-fopenmp
endif

ifneq ($(gpl),0)
	AOBJS+=QSufSort.o bwtgen.o
	override CPPFLAGS+=-DUSE_GPL
endif

ifeq ($(mimalloc),0)
	MALLOC_O=
	override CPPFLAGS+=-DHAVE_KALLOC
endif

ifeq ($(ARCH), x86_64)
	override CFLAGS+=-msse4.2 -mpopcnt
endif

# ===================== B2: compile-optional cp_occ seeding backend =========
# `make B2=1` links b2idx.o (minibwa's seeding leaves re-expressed over
# bwa-mem2/bwa-mem3's cp_occ FM-index) against a pre-built bwa-mem3 libbwa.a.
# b2idx.cpp is C++ (FMI_search is C++) and MUST be compiled with the same ISA
# shim macros / include path bwa-mem3 used to build libbwa.a or the FMI_search
# struct ABI mismatches. Override BWAMEM3_DIR to point at that checkout.
#
# Without B2=1: no b2idx.o, no libbwa.a link, no MB_HAVE_B2 -- native build is
# byte-for-byte the same command line as before this guard existed.
LINK       = $(CC)
LINK_FLAGS = $(CFLAGS)
B2OBJS     =
B3_LDLIBS  =
ifeq ($(B2),1)
	override CPPFLAGS += -DMB_HAVE_B2
	CXX ?= c++
	BWAMEM3_DIR ?= /Users/nhomer/work/git/bwa-mem3/main
	UNAME_S := $(shell uname -s)
	ifneq (,$(filter arm64 aarch64,$(ARCH)))
		B3_ARCH_FLAGS = -DAPPLE_SILICON=1 -DCACHE_LINE_BYTES=128 \
		                -D__SSE__=1 -D__SSE2__=1 -D__SSE3__=1 -D__SSSE3__=1 \
		                -D__SSE4_1__=1 -D__SSE4_2__=1 -I$(BWAMEM3_DIR)/ext/sse2neon
	else
		B3_ARCH_FLAGS =
	endif
	# libbwa.a's FMI_search.o references libsais_build_fm_index(), which in
	# turn needs libsais{,64}_gsa_omp -- present in bwa-mem3's own libsais.o /
	# libsais64.o but NOT baked into libbwa.a itself, so link them explicitly.
	B3_LIBSAIS = $(BWAMEM3_DIR)/ext/libsais/src/libsais.o $(BWAMEM3_DIR)/ext/libsais/src/libsais64.o
	ifeq ($(UNAME_S),Darwin)
		LIBDEFLATE_PREFIX := $(shell brew --prefix libdeflate 2>/dev/null)
		ifneq (,$(wildcard /opt/homebrew/opt/libomp/lib/libomp.dylib))
			B3_OMP = -L/opt/homebrew/opt/libomp/lib -lomp
		else
			B3_OMP = -lomp
		endif
		B3_FRAMEWORKS = -framework Accelerate
		B3_EXTRA = -L$(LIBDEFLATE_PREFIX)/lib -ldeflate
	else
		B3_OMP = -fopenmp
		B3_FRAMEWORKS =
		B3_EXTRA = -ldeflate
		ifneq (,$(wildcard /usr/local/lib64/libdeflate.*))
			B3_EXTRA += -L/usr/local/lib64
		endif
	endif
	# zlib-ng: src/fast_reader.c's plain-gzip decode path is a direct (zng_*)
	# dependency of libbwa.a as of the combined B1+B2 bwa-mem3 -- the older
	# bwa-mem3 checkout the b2idx prototype was validated against predates
	# this, so it isn't in the prototype's Makefile. Use the vendored static
	# archive bwa-mem3 itself builds (ext/zlib-ng/build/libz-ng.a).
	B3_ZLIBNG = $(BWAMEM3_DIR)/ext/zlib-ng/build/libz-ng.a
	B3_LDLIBS = $(BWAMEM3_DIR)/libbwa.a $(BWAMEM3_DIR)/ext/htslib/libhts.a $(B3_LIBSAIS) \
	            $(B3_ZLIBNG) $(B3_OMP) $(B3_FRAMEWORKS) $(B3_EXTRA)
	B2OBJS = b2idx.o
	LINK       = $(CXX)
	LINK_FLAGS =
	# bwa-mem3's libbwa.a pulls in its own libsais (needs the gsa_omp variant
	# above); linking minibwa's libsais.o/libsais64.o too would duplicate
	# symbols, so drop them from AOBJS. (`map` never builds an index, so this
	# only affects the unused `index` subcommand under a B2 build.)
	AOBJS := $(filter-out libsais.o libsais64.o,$(AOBJS))
endif
# ===========================================================================

.SUFFIXES:.c .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

b2idx.o:b2idx.cpp b2idx.h bwt.h
		$(CXX) -c -std=c++14 -O3 -g -Wall $(B3_ARCH_FLAGS) -I$(BWAMEM3_DIR)/src -I. $< -o $@

mimalloc.o:
		$(CC) -c -std=gnu11 -O3 -Wall -Wextra -DNDEBUG -DMI_MALLOC_OVERRIDE -DMI_OSX_INTERPOSE=1 -DMI_OSX_ZONE=1 -Imimalloc mimalloc/static.c -o $@

libminibwa.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

minibwa:libminibwa.a $(MALLOC_O) $(AOBJS) $(B2OBJS) main.o
		$(LINK) $(LINK_FLAGS) $(LDFLAGS) $(MALLOC_O) $(AOBJS) $(B2OBJS) main.o -o $@ -L. -lminibwa $(LIBS) $(B3_LDLIBS)

clean:
		rm -fr *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c *.cpp)

# DO NOT DELETE

QSufSort.o: QSufSort.h
align.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h bseq.h kalloc.h ksw2.h
bseq.o: bseq.h kommon.h kseq.h
bwt.o: kommon.h kalloc.h bwt.h
bwtgen.o: QSufSort.h
cs.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h bseq.h kalloc.h
fastmap.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h bseq.h ketopt.h kseq.h
fastmap.o: kalloc.h
format.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h bseq.h
index.o: libsais.h libsais64.h kommon.h ketopt.h mbpriv.h minibwa.h l2bit.h
index.o: bwt.h bseq.h
kalloc.o: kalloc.h
kommon.o: kommon.h
ksw2_extd2_sse.o: ksw2.h
ksw2_extz2_sse.o: ksw2.h
ksw2_ll_sse.o: ksw2.h
kthread.o: kthread.h
l2bit.o: kommon.h l2bit.h kseq.h
lchain.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h bseq.h kalloc.h ksort.h
libsais.o: libsais.h
libsais64.o: libsais.h libsais64.h
main.o: kommon.h mbpriv.h minibwa.h l2bit.h bwt.h bseq.h ketopt.h
map-algo.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h bseq.h kalloc.h ksort.h
map-main.o: kommon.h mbpriv.h minibwa.h l2bit.h bwt.h bseq.h kalloc.h
map-main.o: kthread.h ketopt.h kseq.h
options.o: minibwa.h
pe.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h bseq.h kalloc.h ksw2.h
seed.o: mbpriv.h minibwa.h l2bit.h bwt.h kommon.h bseq.h kalloc.h ksort.h
