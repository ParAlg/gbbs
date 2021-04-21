# git root directory
ROOTDIR = $(strip $(shell git rev-parse --show-cdup))

include $(ROOTDIR)makefile.variables

ALL= HACSimilarity HACDissimilarity

OTHERFLAGS = -I$(ROOTDIR)external/PAM/include/

include $(ROOTDIR)benchmarks/makefile.benchmarks
