# git root directory
ROOTDIR = $(strip $(shell git rev-parse --show-cdup))

include $(ROOTDIR)makefile.variables

ALL= graph_filter_test

include $(ROOTDIR)benchmarks/makefile.benchmarks

