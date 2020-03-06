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

INCLUDE_DIRS = 

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

ALL= alloc.o get_time.o memory_size.o parallel.o scheduler.o time_operations.o utilities.o

all: $(ALL)

time_operations.o : time_operations.cc get_time.o utilities.o
	$(CC) $(CFLAGS) $(PFLAGS) -c $<

utilities.o : utilities.cc parallel.o
	$(CC) $(CFLAGS) $(PFLAGS) -c $<

alloc.o : alloc.cc memory_size.o parallel.o
	$(CC) $(CFLAGS) $(PFLAGS) -c $<

%.o : %.cc
	$(CC) $(CFLAGS) $(PFLAGS) -c $<

.PHONY : clean

clean :
	rm -f *.o $(ALL)
