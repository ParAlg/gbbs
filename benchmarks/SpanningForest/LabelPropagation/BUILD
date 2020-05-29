load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "SpanningForest",
  hdrs = ["SpanningForest.h"],
  deps = [
  "//benchmarks/SpanningForest:common",
  "//ligra:ligra",
  "//pbbslib:random_shuffle",
  ]
)

gbbs_cc_binary(
  name = "SpanningForest_main",
  srcs = ["SpanningForest.cc"],
  deps = [
  ":SpanningForest",
  "//benchmarks/SpanningForest:common",
  "//benchmarks/SpanningForest/BFSSF:SpanningForest",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
