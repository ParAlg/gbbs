load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "DensestSubgraph",
  hdrs = ["DensestSubgraph.h"],
  deps = [
  "//ligra:edge_map_reduce",
  "//ligra:ligra"
  ]
)

gbbs_cc_binary(
  name = "DensestSubgraph_main",
  srcs = ["DensestSubgraph.cc"],
  deps = [":DensestSubgraph"]
)

package(
  default_visibility = ["//visibility:public"],
)

