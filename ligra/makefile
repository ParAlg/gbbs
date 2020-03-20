include ../makefile.variables

INCLUDE_DIRS = -I../

OBJDIR = ../bin/ligra/

ALL_PRE = bridge edge_map_blocked graph graph_io io ligra parse_command_line undirected_edge union_find
ALL= $(addprefix $(OBJDIR), $(addsuffix .o, $(ALL_PRE))) $(addprefix $(OBJDIR), $(addsuffix .a, $(ALL_PRE)))

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

$(OBJDIR)%.a : $(OBJDIR)%.o
	libtool -static -s -o $@ $<

.PHONY : clean

clean :
	rm -f *.o $(ALL)
