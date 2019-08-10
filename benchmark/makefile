ifeq (, $(shell which jemalloc-config))
JEMALLOC =
else
JEMALLOCLD = $(shell jemalloc-config --libdir)
JEMALLOC = -L$(JEMALLOCLD) -ljemalloc
endif

# Always compile with LONG (note that timings on small graphs may be a bit
# faster w/o this flag).
INTT = -DLONG

ifdef EDGELONG
INTE = -DEDGELONG
endif

#INCLUDE_DIRS = -I/usr0/home/ldhulipa/
INCLUDE_DIRS = -I../

# abm flag on intel indicates lzcnt
# see https://en.wikipedia.org/wiki/Bit_Manipulation_Instruction_Sets
OPT = -O3 -g -mlzcnt -DNDEBUG

ifdef NVM
NVMFLAGS = -lpmem -lnuma -DNVM
endif

CFLAGS = $(INCLUDE_DIRS) -I../src -mcx16 -ldl -std=c++17 -march=native -Wall $(OPT) $(INTT) $(INTE) -DAMORTIZEDPD $(CONCEPTS) -DUSEMALLOC $(NVMFLAGS)

OMPFLAGS = -DOPENMP -fopenmp
CILKFLAGS = -DCILK -fcilkplus
HGFLAGS = -DHOMEGROWN -pthread

ifdef CLANG
CC = clang++
PFLAGS = $(CILKFLAGS)
else ifdef CILK
CC = g++
PFLAGS = $(CILKFLAGS)
else ifdef OPENMP
CC = g++
PFLAGS = $(OMPFLAGS)
else ifdef HOMEGROWN
CC = g++
PFLAGS = $(HGFLAGS)
else ifdef SERIAL
CC = g++
PFLAGS =
else # default is homegrown
CC = g++
PFLAGS = $(HGFLAGS)
endif

ALL= BC BellmanFord BFS Biconnectivity CC Coloring DensestSubgraph KCore LDD MaximalMatching MIS MST PageRank SetCover Statistics Spanner SpanningForest Triangle wBFS WidestPath PrimitiveBench

all: $(ALL)

% : %.C
	$(CC) $(CFLAGS) $(PFLAGS) -o $@ $<

.PHONY : clean

clean :
	rm -f *.o $(ALL)
