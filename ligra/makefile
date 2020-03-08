include ../makefile.variables

INCLUDE_DIRS = -I../

ALL= bridge.o edge_map_blocked.o graph.o graph_io.o io.o ligra.o parse_command_line.o undirected_edge.o union_find.o

all: $(ALL)

edge_map_blocked.o : edge_map_blocked.cc bridge.o graph.o
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $<

graph_io.o : graph_io.cc graph.o io.o
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $<

graph.o : graph.cc bridge.o
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $<

ligra.o : ligra.cc bridge.o graph.o graph_io.o parse_command_line.o
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $<

%.o : %.cc
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $<

.PHONY : clean

clean :
	rm -f *.o $(ALL)
