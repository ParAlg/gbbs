include ../../makefile.variables

INCLUDE_DIRS = -I../../

OBJDIR = ../../bin/ligra/encodings/

ALL= $(addprefix $(OBJDIR), byte.o byte_pd_amortized.o)

all: $(ALL)

$(OBJDIR)byte.o : byte.cc
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)byte_pd_amortized.o : byte_pd_amortized.cc
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $< -o $@

.PHONY : clean

clean :
	rm -f *.o $(ALL)
