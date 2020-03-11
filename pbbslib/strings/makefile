include ../../makefile.variables

ALL= ../../bin/pbbslib/strings/string_basics.o

all: $(ALL)

../../bin/pbbslib/strings/string_basics.o : string_basics.cc
	$(CC) $(CFLAGS) $(PFLAGS) -c $< -o $@

.PHONY : clean

clean :
	rm -f *.o $(ALL)
