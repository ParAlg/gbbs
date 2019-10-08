cc_library(
  name = "DensestSubgraph",
  hdrs = ["DensestSubgraph.h"],
  deps = [
  "//ligra:edge_map_reduce",
  "//ligra:ligra"
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

