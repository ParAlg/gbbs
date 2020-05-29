load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "contract_sf",
  hdrs = ["contract_sf.h"],
  deps = [
  "//ligra:graph",
  "//ligra:macros",
  "//ligra/pbbslib:sparse_table",
  "//ligra/pbbslib:dyn_arr",
  "//pbbslib:sequence_ops"
  ]
)

gbbs_cc_library(
  name = "SpanningForest",
  hdrs = ["SpanningForest.h"],
  deps = [
  ":contract_sf",
  "//benchmarks/LowDiameterDecomposition/MPX13:LowDiameterDecomposition",
  "//ligra:ligra",
  "//ligra/pbbslib:sparse_table"
  ]
)

gbbs_cc_binary(
  name = "SpanningForest_main",
  srcs = ["SpanningForest.cc"],
  deps = [
  ":SpanningForest",
  "//benchmarks/SpanningForest:common",
  "//benchmarks/SpanningForest/BFSSF:SpanningForest"
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
