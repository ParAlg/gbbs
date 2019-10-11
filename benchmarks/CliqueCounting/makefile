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
INCLUDE_DIRS = -I../../ -I../../external/simdinter/include

OPT = -O3 -g

CFLAGS = $(INCLUDE_DIRS) -mcx16 -ldl -std=c++17 -march=native $(OPT) $(INTT) $(INTE) -DAMORTIZEDPD $(CONCEPTS) -DUSEMALLOC
#SIMDFLAGS = -I../simdinter -I../simdinter/include -fpic -mavx -mcx16 -ldl -std=c++17 -march=native $(OPT) $(INTT) $(INTE) -DAMORTIZEDPD $(CONCEPTS) -DUSEMALLOC

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

ALL= Clique

all: $(ALL)

#intersection.o: ../simdinter/include/intersection.h ../simdinter/src/intersection.cpp
#	$(CC) $(SIMDFLAGS) -c ../simdinter/src/intersection.cpp -Iinclude

#OBJECTS = intersection.o

% : %.cc simdinter graphsetinter
	$(CC) $(CFLAGS) $(PFLAGS) -o $@ $< ../../external/simdinter/intersection.o ../../external/graphsetinter/src/set_operation.o

.PHONY : clean simdinter

simdinter:
	make intersection.o -C ../../external/simdinter

.PHONY : clean graphsetinter

graphsetinter:
	make set_operation.o -C ../../external/graphsetinter/src

clean :
	rm -f *.o $(ALL)
