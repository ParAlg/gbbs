include ../makefile.variables

INCLUDE_DIRS = -I../ -I../parlaylib/include/

OBJDIR = ../bin/gbbs/

$(info $$var is [$(ROOTDIR)])
$(info VAR is $(OBJDIR))
$(info VAR is $(ROOTDIR))
$(info VAR is $(CURDIR))

ALL_PRE = benchmark bridge edge_map_blocked get_time graph_io io parse_command_line undirected_edge union_find
ALL= $(addprefix $(OBJDIR), $(addsuffix .o, $(ALL_PRE))) $(addprefix $(OBJDIR), $(addsuffix .a, $(ALL_PRE)))

all: $(ALL)

$(OBJDIR)edge_map_blocked.o : edge_map_blocked.cc $(addprefix $(OBJDIR), bridge.o)
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)graph_io.o : graph_io.cc $(addprefix $(OBJDIR), io.o)
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)%.o : %.cc
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)%.a : $(OBJDIR)%.o
	ar -rcs $@ $<

.PHONY : clean

clean :
	rm -f *.o $(ALL)
