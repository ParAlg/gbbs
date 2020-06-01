cc_library(
  name = "DensestSubgraph",
  hdrs = ["DensestSubgraph.h"],
  deps = [
  "//gbbs:edge_map_reduce",
  "//gbbs:gbbs"
  ]
)

cc_binary(
  name = "DensestSubgraph_main",
  srcs = ["DensestSubgraph.cc"],
  deps = [
    ":DensestSubgraph",
    "//sage:sage",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
