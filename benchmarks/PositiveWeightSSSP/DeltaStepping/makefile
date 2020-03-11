# git root directory
ROOTDIR = $(strip $(shell git rev-parse --show-cdup))

include $(ROOTDIR)makefile.variables

ALL= DeltaStepping

include $(ROOTDIR)benchmarks/makefile.benchmarks

