include ../../makefile.variables

ALL= ../../bin/pbbslib/strings/string_basics.o ../../bin/pbbslib/strings/string_basics.a

all: $(ALL)

../../bin/pbbslib/strings/string_basics.o : string_basics.cc
	$(CC) $(CFLAGS) $(PFLAGS) -c $< -o $@

../../bin/pbbslib/strings/string_basics.a : ../../bin/pbbslib/strings/string_basics.o
	libtool -static -s -o $@ $<

.PHONY : clean

clean :
	rm -f *.o $(ALL)
