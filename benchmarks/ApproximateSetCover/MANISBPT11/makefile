# git root directory
ROOTDIR = $(strip $(shell git rev-parse --show-cdup))

include $(ROOTDIR)makefile.variables

ALL= ApproximateSetCover

include $(ROOTDIR)benchmarks/makefile.benchmarks
