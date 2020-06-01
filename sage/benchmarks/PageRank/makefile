# git root directory
ROOTDIR = $(strip $(shell git rev-parse --show-cdup))

include $(ROOTDIR)makefile.variables

ALL= PageRank

include $(ROOTDIR)benchmarks/makefile.benchmarks

