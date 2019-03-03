CFLAGS = -I../ -mcx16 -O3 -std=c++17 -Wall

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
else # default is cilk
CC = g++
PFLAGS = $(CILKFLAGS)
endif

EXAMPLES = mcss

# object files
% : %.cpp
	$(CC) $(CFLAGS) $(PFLAGS) $@.cpp -o $@

clean:
	rm -f $(EXAMPLES)
