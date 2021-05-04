# git root directory
ROOTDIR = $(strip $(shell git rev-parse --show-cdup))

include $(ROOTDIR)makefile.variables

ALL= ParallelHAC

OTHERFLAGS = -I$(ROOTDIR)external/PAM/include/

include $(ROOTDIR)benchmarks/makefile.benchmarks
