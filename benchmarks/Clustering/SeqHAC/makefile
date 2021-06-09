# git root directory
ROOTDIR = $(strip $(shell git rev-parse --show-cdup))

include $(ROOTDIR)makefile.variables

#ALL= HACSimilarity HACDissimilarity
ALL= HACSimilarity

OTHERFLAGS = -I$(ROOTDIR)external/PAM/include/ -DPARLAY_SEQUENTIAL

include $(ROOTDIR)benchmarks/makefile.benchmarks
