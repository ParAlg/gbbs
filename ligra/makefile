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

INCLUDE_DIRS = -I../

OPT = -O3 -g -DNDEBUG

CFLAGS = $(INCLUDE_DIRS) -mcx16 -std=c++17 -march=native -Wall $(OPT) $(INTT) $(INTE) -DAMORTIZEDPD $(CONCEPTS) -DUSEMALLOC -DNDEBUG

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

ALL= bridge.o graph.o graph_io.o io.o ligra.o parse_command_line.o undirected_edge.o union_find.o

all: $(ALL)

.PHONY : clean pbbslib

pbbslib:
	make -C ../pbbslib/

graph_io.o : graph_io.cc graph.o io.o pbbslib
	$(CC) $(CFLAGS) $(PFLAGS) -c $<

graph.o : graph.cc bridge.o pbbslib
	$(CC) $(CFLAGS) $(PFLAGS) -c $<

ligra.o : ligra.cc bridge.o graph.o graph_io.o parse_command_line.o pbbslib
	$(CC) $(CFLAGS) $(PFLAGS) -c $<

%.o : %.cc pbbslib
	$(CC) $(CFLAGS) $(PFLAGS) -c $<

.PHONY : clean

clean :
	rm -f *.o $(ALL)
