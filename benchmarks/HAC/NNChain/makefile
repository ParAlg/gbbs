# git root directory
ROOTDIR = $(strip $(shell git rev-parse --show-cdup))

include $(ROOTDIR)makefile.variables

ALL= HAC

#OTHER_OBJS = \
#  $(wildcard $(ROOTDIR)bin/benchmarks/Connectivity/*.o)
#OTHER = \
#  connectivity_objs

## For `benchmarks/Connectivity/UnionFind` dependency
#connectivity_objs:
#	make -C $(ROOTDIR)benchmarks/Connectivity/

include $(ROOTDIR)benchmarks/makefile.benchmarks
