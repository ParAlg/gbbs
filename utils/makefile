ifdef LONG
INTT = -DLONG
endif

ifdef EDGELONG
INTE = -DEDGELONG
endif

PD = -DAMORTIZEDPD

PCC = g++
PCFLAGS = -std=c++14 -mcx16 -fcilkplus -lcilkrts -I../ -I../src -O2 -DCILK $(INTT) $(INTE) $(CODE) $(PD) $(MEM)
PLFLAGS = -fcilkplus -lcilkrts

ALL= add_weights converter

all: $(ALL)

% : %.C
	$(PCC) $(PCFLAGS) -o $@ $<

.PHONY : clean

clean :
	rm -f *.o $(ALL)

cleansrc :
	rm -f *.o $(ALL)
	rm $(COMMON)
	rm -rf lib/*.h
