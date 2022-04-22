include ../../makefile.variables

INCLUDE_DIRS = -I../../ -I../../external/parlaylib/include/

OBJDIR = ../../bin/benchmarks/Connectivity/

ALL= $(addprefix $(OBJDIR), common.o connectit.o)

all: $(ALL)

$(OBJDIR)common.o : common.cc
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)connectit.o : connectit.cc
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $< -o $@

.PHONY : clean

clean :
	rm -f *.o $(ALL)
