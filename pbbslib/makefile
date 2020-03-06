include ../makefile.variables

ALL= alloc.o get_time.o memory_size.o parallel.o scheduler.o time_operations.o utilities.o

all: $(ALL)

time_operations.o : time_operations.cc get_time.o utilities.o
	$(CC) $(CFLAGS) $(PFLAGS) -c $<

utilities.o : utilities.cc parallel.o
	$(CC) $(CFLAGS) $(PFLAGS) -c $<

alloc.o : alloc.cc memory_size.o parallel.o
	$(CC) $(CFLAGS) $(PFLAGS) -c $<

%.o : %.cc
	$(CC) $(CFLAGS) $(PFLAGS) -c $<

.PHONY : clean

clean :
	rm -f *.o $(ALL)
