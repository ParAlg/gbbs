cc_library(
  name = "DensestSubgraph",
  hdrs = ["DensestSubgraph.h"],
  deps = [
  "//gbbs:edge_map_reduce",
  "//gbbs:gbbs",
  "//benchmarks/KCore/JulienneDBS17:KCore"
  ]
)

cc_binary(
  name = "DensestSubgraph_main",
  srcs = ["DensestSubgraph.cc"],
  deps = [":DensestSubgraph"]
)

package(
  default_visibility = ["//visibility:public"],
)

