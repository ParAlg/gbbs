ifdef LONG
INTT = -DLONG
endif

ifdef EDGELONG
INTE = -DEDGELONG
endif

OPT = -O2

INCLUDE_DIRS = -I/usr0/home/ldhulipa/
#INCLUDE_DIRS = -I../

PCC = g++
PCFLAGS = -std=c++14 -fcilkplus -lcilkrts -I../src $(INCLUDE_DIRS) -mcx16 -march=native $(OPT) -DHOMEGROWN -pthread -DAMORTIZEDPD $(INTT) $(INTE)
#PCFLAGS = -std=c++14 -fcilkplus -lcilkrts -I../src $(INCLUDE_DIRS) -mcx16 -march=native $(OPT) -DCILK -DAMORTIZEDPD $(INTT) $(INTE)

ALL= BC BellmanFord BFS Biconnectivity CC Coloring DensestSubgraph KCore KTruss LDD MaximalMatching MIS MST SCC SetCover Triangle wBFS

all: $(ALL)

% : %.C
	$(PCC) $(PCFLAGS) -o $@ $<

.PHONY : clean

clean :
	rm -f *.o $(ALL)
