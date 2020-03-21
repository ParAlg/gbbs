include ../makefile.variables

OBJDIR = ../bin/pbbslib/

ALL_PRE = alloc get_time memory_size parallel scheduler time_operations utilities
ALL= $(addprefix $(OBJDIR), $(addsuffix .o, $(ALL_PRE))) $(addprefix $(OBJDIR), $(addsuffix .a, $(ALL_PRE)))

all: $(ALL)

$(OBJDIR)time_operations.o : time_operations.cc $(addprefix $(OBJDIR), get_time.o utilities.o)
	$(CC) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)utilities.o : utilities.cc $(addprefix $(OBJDIR), parallel.o)
	$(CC) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)alloc.o : alloc.cc $(addprefix $(OBJDIR), memory_size.o parallel.o)
	$(CC) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)%.o : %.cc
	$(CC) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)%.a : $(OBJDIR)%.o
	ar -rcs $@ $<

.PHONY : clean

clean :
	rm -f *.o $(ALL)
