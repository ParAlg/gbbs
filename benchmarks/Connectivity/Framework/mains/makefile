ifeq (, $(shell which jemalloc-config))
JEMALLOC =
else
JEMALLOCLD = $(shell jemalloc-config --libdir)
JEMALLOC = -L$(JEMALLOCLD) -ljemalloc
endif

# Always compile with LONG (note that timings on small graphs with less than
# 2^{32} edges may be slightly faster without this flag)
INTT = -DLONG

ifdef EDGELONG
INTE = -DEDGELONG
endif

FORMAT = -DCSR # (oneof {DCSR, DCOO})

INCLUDE_DIRS = -I../../../../

#OPT = -O3 -g
OPT = -O3 -g -DNDEBUG
#OPT = -O0 -g

#INTEL_PCM_HOME = /usr0/home/ldhulipa/pcm/

#ifdef INTEL_PCM_HOME
#	INCLUDE_DIRS+=-I$(INTEL_PCM_HOME)
#	EXT_OBJECTS=$(INTEL_PCM_HOME)/msr.o $(INTEL_PCM_HOME)/cpucounters.o $(INTEL_PCM_HOME)/pci.o $(INTEL_PCM_HOME)/client_bw.o $(INTEL_PCM_HOME)/utils.o $(INTEL_PCM_HOME)/mmio.o
#	PCM = -DUSE_PCM_LIB
#endif

CFLAGS = $(INCLUDE_DIRS) -mcx16 -ldl -std=c++17 -march=native $(OPT) $(INTT) $(INTE) -DAMORTIZEDPD $(CONCEPTS) -DUSEMALLOC $(FORMAT) $(PCM) $(EXT_OBJECTS)

OMPFLAGS = -DOPENMP -fopenmp
CILKFLAGS = -DCILK -fcilkplus
HGFLAGS = -DHOMEGROWN -pthread

ifdef CLANG
CC = clang++
PFLAGS = $(CILKFLAGS)
else ifdef CILK
CC = g++-7
PFLAGS = $(CILKFLAGS)
else ifdef OPENMP
CC = g++-7
PFLAGS = $(OMPFLAGS)
else ifdef HOMEGROWN
CC = g++-7
PFLAGS = $(HGFLAGS)
else ifdef SERIAL
CC = g++-7
PFLAGS =
else # default is homegrown
CC = g++-7
PFLAGS = $(HGFLAGS)
endif

ALL=unite_nosample unite_nd_nosample unite_early_nosample unite_rem_lock_nosample unite_rem_cas_nosample unite_kout unite_nd_kout unite_early_kout unite_rem_lock_kout unite_rem_cas_kout unite_bfs unite_nd_bfs unite_early_bfs unite_rem_lock_bfs unite_rem_cas_bfs unite_ldd unite_nd_ldd unite_early_ldd unite_rem_lock_ldd unite_rem_cas_ldd jayanti_nosample jayanti_kout jayanti_bfs jayanti_ldd label_propagation liutarjan_nosample shiloach_vishkin bfscc gbbscc

all: $(ALL)

% : %.cc
	$(CC) $(CFLAGS) $(PFLAGS) -o $@ $<

.PHONY : clean

clean :
	rm -f *.o $(ALL)
