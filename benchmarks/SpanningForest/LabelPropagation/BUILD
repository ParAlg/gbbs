licenses(["notice"])

cc_library(
  name = "SpanningForest",
  hdrs = ["SpanningForest.h"],
  deps = [
  "//gbbs:gbbs",
  "//benchmarks/Connectivity:common",
  "//benchmarks/SpanningForest:common",
  ]
)

cc_binary(
  name = "SpanningForest_main",
  srcs = ["SpanningForest.cc"],
  deps = [
  ":SpanningForest",
  "//gbbs:gbbs",
  "//benchmarks/SpanningForest:common",
  "//benchmarks/SpanningForest/BFSSF:SpanningForest",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
