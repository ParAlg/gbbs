# git root directory
ROOTDIR = $(strip $(shell git rev-parse --show-cdup))

include $(ROOTDIR)makefile.variables

ALL= Connectivity
OTHER = connectivity_objs

OTHER_OBJS = $(wildcard $(ROOTDIR)bin/benchmarks/Connectivity/*.o)

connectivity_objs :
	make -C $(ROOTDIR)benchmarks/Connectivity/

include $(ROOTDIR)benchmarks/makefile.benchmarks
