ifdef LONG
INTT = -DLONG
endif

ifdef EDGELONG
INTE = -DEDGELONG
endif

OPT = -O3

#INCLUDE_DIRS = -I/usr0/home/ldhulipa/
INCLUDE_DIRS = -I../

PCC = g++
PCFLAGS = -std=c++14 -fcilkplus -lcilkrts -I../src $(INCLUDE_DIRS) -mcx16 -march=native $(OPT) -DCILK -DAMORTIZEDPD $(INTT) $(INTE)

ALL= BC BellmanFord BFS Biconnectivity CC Coloring KCore LDD MaximalMatching MIS MST SCC SetCover Triangle wBFS

all: $(ALL)

% : %.C
	$(PCC) $(PCFLAGS) -o $@ $<

.PHONY : clean

clean :
	rm -f *.o $(ALL)
