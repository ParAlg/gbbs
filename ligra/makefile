include ../makefile.variables

ALL= bridge.o graph.o graph_io.o io.o ligra.o parse_command_line.o undirected_edge.o union_find.o

all: $(ALL)

.PHONY : clean pbbslib

pbbslib:
	make -C ../pbbslib/

graph_io.o : graph_io.cc graph.o io.o pbbslib
	$(CC) $(CFLAGS) $(PFLAGS) -c $<

graph.o : graph.cc bridge.o pbbslib
	$(CC) $(CFLAGS) $(PFLAGS) -c $<

ligra.o : ligra.cc bridge.o graph.o graph_io.o parse_command_line.o pbbslib
	$(CC) $(CFLAGS) $(PFLAGS) -c $<

%.o : %.cc pbbslib
	$(CC) $(CFLAGS) $(PFLAGS) -c $<

.PHONY : clean

clean :
	rm -f *.o $(ALL)
