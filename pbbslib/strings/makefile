ifeq (, $(shell which jemalloc-config))
JEMALLOC =
else
JEMALLOCLD = $(shell jemalloc-config --libdir)
JEMALLOC = -L$(JEMALLOCLD) -ljemalloc 
endif

CONCEPTS = -fconcepts -DCONCEPTS
CFLAGS = -mcx16 -O3 -ldl -std=c++17 -march=native -Wall 

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

INCLUDE = -I../

test_suffix_tree: suffix_tree.h suffix_array.h lcp.h test_suffix_tree.cpp
	$(CC) $(CFLAGS) $(PFLAGS) $(INCLUDE) test_suffix_tree.cpp -o test_suffix_tree $(JEMALLOC)
