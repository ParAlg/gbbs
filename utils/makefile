# Builds binaries in utils/.

# git root directory
ROOTDIR = $(strip $(shell git rev-parse --show-cdup))

include $(ROOTDIR)makefile.variables

ALL = \
	compressor \
	converter \
	random_reorder \
	to_edge_list \
	snap_converter

# These aren't benchmarks, but they can be built in the same fashion.
include $(ROOTDIR)benchmarks/makefile.benchmarks
