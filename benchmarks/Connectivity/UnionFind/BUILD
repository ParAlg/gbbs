load("//internal_tools:build_defs.bzl", "gbbs_cc_library")

gbbs_cc_library(
  name = "union_find_rules",
  hdrs = ["union_find_rules.h"],
  deps = [
  ":jayanti",
  "//ligra:bridge",
  "//pbbslib:seq",
  "//benchmarks/Connectivity:common",
  ]
)

gbbs_cc_library(
  name = "jayanti",
  hdrs = ["jayanti.h"],
  deps = [
  "//ligra:ligra",
  "//pbbslib:seq",
  "//benchmarks/Connectivity:common",
  ]
)

gbbs_cc_library(
  name = "Connectivity",
  hdrs = ["Connectivity.h"],
  deps = [
  ":union_find_rules",
  "//ligra:bridge",
  "//ligra:ligra",
  "//ligra/pbbslib:sparse_table"
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
