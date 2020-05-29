load("//internal_tools:build_defs.bzl", "gbbs_cc_library")

gbbs_cc_library(
  name = "liu_tarjan_rules",
  hdrs = ["liu_tarjan_rules.h"],
  deps = [
  "//benchmarks/Connectivity:common",
  ]
)

gbbs_cc_library(
  name = "Connectivity",
  hdrs = ["Connectivity.h"],
  deps = [
  ":liu_tarjan_rules",
  "//ligra:ligra",
  "//ligra:bridge",
  "//pbbslib:random",
  "//benchmarks/Connectivity:common",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
