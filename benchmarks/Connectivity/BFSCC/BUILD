load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "Connectivity",
  hdrs = ["Connectivity.h"],
  deps = [
  "//ligra:ligra",
  "//benchmarks/Connectivity:common",
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
