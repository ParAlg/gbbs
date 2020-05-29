load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "Spanner",
  hdrs = ["Spanner.h"],
  deps = [
  "//ligra:ligra",
  "//ligra/pbbslib:sparse_table",
  "//ligra/pbbslib:dyn_arr",
  "//benchmarks/LowDiameterDecomposition/MPX13:LowDiameterDecomposition",
  ]
)

gbbs_cc_binary(
  name = "Spanner_main",
  srcs = ["Spanner.cc"],
  deps = [":Spanner"]
)

package(
  default_visibility = ["//visibility:public"],
)
