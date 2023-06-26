CC=			gcc
CFLAGS=		-g -Wall -Wno-unused-function -O2
WRAP_MALLOC=
AR=			ar
DFLAGS=		-DHAVE_PTHREAD $(WRAP_MALLOC)
LOBJS=		bwalib/utils.o \
			bwalib/ksw.o \
			bwalib/bwa.o \
			cstl/kthread.o \
			cstl/kstring.o \
			FM_index/bwt.o \
			FM_index/bntseq.o \
			FM_index/QSufSort.o \
			FM_index/bwt_gen.o \
			FM_index/rope.o \
			FM_index/rle.o \
			FM_index/is.o
AOBJS=		bwalib/bwashm.o \
			bwalib/kopen.o
PROG=		bwamem_seeding
INCLUDES=
LIBS=		-lm -lz -lpthread

ifeq ($(shell uname -s),Linux)
	LIBS += -lrt
endif

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

$(PROG):libbwa.a $(AOBJS) seeding/bwa_seeding.o
		$(CC) $(CFLAGS) $(DFLAGS) $(AOBJS) seeding/bwa_seeding.o -o $@ -L. -lbwa $(LIBS)

libbwa.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

clean:
		rm -f FM_index/*.o cstl/*.o bwalib/*.o libbwa.a $(PROG)
