# Builds binaries in utils/generators.
#
# git root directory
ROOTDIR = $(strip $(shell git rev-parse --show-cdup))

include $(ROOTDIR)makefile.variables

ALL = \
	RMAT \
	BarabasiAlbert

# These aren't benchmarks, but they can be built in the same fashion.
include $(ROOTDIR)benchmarks/makefile.benchmarks
