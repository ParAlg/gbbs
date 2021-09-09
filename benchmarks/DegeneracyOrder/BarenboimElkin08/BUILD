licenses(["notice"])

cc_library(
  name = "DegeneracyOrder",
  hdrs = ["DegeneracyOrder.h"],
  deps = [
  "//gbbs:bucket",
  "//gbbs:gbbs",
  "//gbbs/helpers:dyn_arr",
  "//benchmarks/ApproximateDensestSubgraph/GreedyCharikar:DensestSubgraph",
  "//benchmarks/ApproximateDensestSubgraph/ApproxPeelingBKV12:DensestSubgraph",
  ]
)

cc_binary(
  name = "DegeneracyOrder_main",
  srcs = ["DegeneracyOrder.cc"],
  deps = [":DegeneracyOrder"]
)

package(
  default_visibility = ["//visibility:public"],
)
