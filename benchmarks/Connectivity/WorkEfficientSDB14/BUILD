load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "contract",
  hdrs = ["contract.h"],
  deps = [
  "//ligra:graph",
  "//ligra:macros",
  "//ligra/pbbslib:sparse_table",
  "//pbbslib:sequence_ops"
  ]
)

gbbs_cc_library(
  name = "Connectivity",
  hdrs = ["Connectivity.h"],
  deps = [
  ":contract",
  "//benchmarks/LowDiameterDecomposition/MPX13:LowDiameterDecomposition",
  "//benchmarks/Connectivity:common",
  "//ligra:ligra",
  "//ligra/pbbslib:sparse_table"
  ]
)

gbbs_cc_binary(
  name = "Connectivity_main",
  srcs = ["Connectivity.cc"],
  deps = [":Connectivity"]
)

package(
  default_visibility = ["//visibility:public"],
)
