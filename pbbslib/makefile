include ../makefile.variables

OBJDIR = ../bin/pbbslib/

ALL= $(addprefix $(OBJDIR), alloc.o get_time.o memory_size.o parallel.o scheduler.o time_operations.o utilities.o)

all: $(ALL)

$(OBJDIR)time_operations.o : time_operations.cc $(addprefix $(OBJDIR), get_time.o utilities.o)
	$(CC) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)utilities.o : utilities.cc $(addprefix $(OBJDIR), parallel.o)
	$(CC) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)alloc.o : alloc.cc $(addprefix $(OBJDIR), memory_size.o parallel.o)
	$(CC) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)%.o : %.cc
	$(CC) $(CFLAGS) $(PFLAGS) -c $< -o $@

.PHONY : clean

clean :
	rm -f *.o $(ALL)
