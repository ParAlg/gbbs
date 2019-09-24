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
INCLUDE_DIRS = -I../ligra -I../

OPT = -O3 -g

CFLAGS = $(INCLUDE_DIRS) -mcx16 -ldl -std=c++17 -march=native -Wall $(OPT) $(INTT) $(INTE) -DAMORTIZEDPD $(CONCEPTS) -DUSEMALLOC
SIMDFLAGS = -I../simdinter -I../simdinter/include -fpic -mavx -mcx16 -ldl -std=c++17 -march=native -Wall $(OPT) $(INTT) $(INTE) -DAMORTIZEDPD $(CONCEPTS) -DUSEMALLOC

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

ALL= AppKCore

all: $(ALL)

#intersection.o: ../simdinter/include/intersection.h ../simdinter/src/intersection.cpp
#	$(CC) $(SIMDFLAGS) -c ../simdinter/src/intersection.cpp -Iinclude

#OBJECTS = intersection.o

% : %.C intersection
	$(CC) $(CFLAGS) $(PFLAGS) -o $@ $< ../simdinter/intersection.o

.PHONY : clean intersection

intersection:
	make intersection.o -C ../simdinter

clean :
	rm -f *.o $(ALL)
