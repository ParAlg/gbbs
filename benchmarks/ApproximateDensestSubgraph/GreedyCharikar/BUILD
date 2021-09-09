licenses(["notice"])

cc_library(
  name = "DensestSubgraph",
  hdrs = ["DensestSubgraph.h"],
  deps = [
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

