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

INCLUDE_DIRS = -I../../

OPT = -O3 -g -DNDEBUG

GOOGLE_TEST = -lgtest -lgtest_main

CFLAGS = $(INCLUDE_DIRS) -mcx16 -ldl -std=c++17 -Wall -I -march=native $(OPT) $(INTT) $(INTE) -DAMORTIZEDPD $(CONCEPTS) -DUSEMALLOC

OMPFLAGS = -DOPENMP -fopenmp
CILKFLAGS = -DCILK -fcilkplus
HGFLAGS = -DHOMEGROWN -pthread

# compile serially for tests
CC = g++-7
PFLAGS =

ALL=graph_test

all: $(ALL)

% : %.cc
	$(CC) $(CFLAGS) $(PFLAGS) -o $@ $< $(GOOGLE_TEST)

.PHONY : clean

clean :
	rm -f *.o $(ALL)
