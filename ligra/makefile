include ../makefile.variables

INCLUDE_DIRS = -I../

OBJDIR = ../bin/ligra/

ALL= $(addprefix $(OBJDIR), bridge.o edge_map_blocked.o graph.o graph_io.o io.o ligra.o parse_command_line.o undirected_edge.o union_find.o)

all: $(ALL)

$(OBJDIR)edge_map_blocked.o : edge_map_blocked.cc $(addprefix $(OBJDIR), bridge.o graph.o)
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)graph_io.o : graph_io.cc $(addprefix $(OBJDIR), graph.o io.o)
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)graph.o : graph.cc $(addprefix $(OBJDIR), bridge.o)
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)ligra.o : ligra.cc $(addprefix $(OBJDIR), bridge.o graph.o graph_io.o parse_command_line.o)
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)%.o : %.cc
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $< -o $@

.PHONY : clean

clean :
	rm -f *.o $(ALL)
