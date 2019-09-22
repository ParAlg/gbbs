cc_library(
  name = "DensestSubgraph",
  hdrs = ["DensestSubgraph.h"],
  deps = [
  "//ligra:edge_map_reduce",
  "//ligra:ligra",
  "//benchmarks/kcore/julienne/KCore.h"
  ]
)

cc_binary(
  name = "DensestSubgraph_main",
  srcs = ["DensestSubgraph.cpp"],
  deps = [":DensestSubgraph"]
)
