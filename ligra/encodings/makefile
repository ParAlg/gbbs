include ../../makefile.variables

INCLUDE_DIRS = -I../../

ALL= byte.o byte_pd_amortized.o

all: $(ALL)

byte.o : byte.cc
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $<

byte_pd_amortized.o : byte_pd_amortized.cc
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $<

.PHONY : clean

clean :
	rm -f *.o $(ALL)
