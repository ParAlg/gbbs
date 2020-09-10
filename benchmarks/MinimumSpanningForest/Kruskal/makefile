# git root directory
ROOTDIR = $(strip $(shell git rev-parse --show-cdup))

include $(ROOTDIR)makefile.variables

ALL = MinimumSpanningForest

MinimumSpanningForest_OBJDIR = $(ROOTDIR)bin/benchmarks/MinimumSpanningForest/Kruskal/
OTHER_LOCAL =
OTHER_OBJS = \
  $(wildcard $(ROOTDIR)bin/benchmarks/Connectivity/*.o) \
  $(addprefix $(MinimumSpanningForest_OBJDIR), $(OTHER_LOCAL))
OTHER = \
  connectivity_objs \
  $(addprefix $(MinimumSpanningForest_OBJDIR), $(OTHER_LOCAL))

# For `benchmarks/Connectivity/UnionFind` dependency
connectivity_objs:
	make -C $(ROOTDIR)benchmarks/Connectivity/

$(MinimumSpanningForest_OBJDIR):
	mkdir -p $@

$(MinimumSpanningForest_OBJDIR)%.o: %.cc %.h | $(MinimumSpanningForest_OBJDIR)
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $< -o $@

include $(ROOTDIR)benchmarks/makefile.benchmarks
