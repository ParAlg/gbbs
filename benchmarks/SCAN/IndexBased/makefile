# git root directory
ROOTDIR = $(strip $(shell git rev-parse --show-cdup))

include $(ROOTDIR)makefile.variables

ALL = run_scan

SCAN_OBJDIR = $(ROOTDIR)bin/benchmarks/SCAN/IndexBased/
OTHER_LOCAL = scan.o scan_helpers.o similarity_measure.o
OTHER_OBJS = \
  $(wildcard $(ROOTDIR)bin/benchmarks/Connectivity/*.o) \
  $(addprefix $(SCAN_OBJDIR), $(OTHER_LOCAL))
OTHER = \
  connectivity_objs \
  $(addprefix $(SCAN_OBJDIR), $(OTHER_LOCAL))

# For `benchmarks/Connectivity/UnionFind` dependency
connectivity_objs:
	make -C $(ROOTDIR)benchmarks/Connectivity/

$(SCAN_OBJDIR):
	mkdir -p $@

$(SCAN_OBJDIR)%.o: %.cc %.h | $(SCAN_OBJDIR)
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $< -o $@

include $(ROOTDIR)benchmarks/makefile.benchmarks
