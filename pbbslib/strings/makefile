include ../../makefile.variables

ALL= string_basics.o

all: $(ALL)

string_basics.o : string_basics.cc
	$(CC) $(CFLAGS) $(PFLAGS) -c $<

.PHONY : clean

clean :
	rm -f *.o $(ALL)
