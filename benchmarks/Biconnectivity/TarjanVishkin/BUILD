load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "Biconnectivity",
  hdrs = ["Biconnectivity.h"],
  deps = [
  "//benchmarks/Connectivity/WorkEfficientSDB14:Connectivity",
  "//ligra:bridge",
  "//ligra:ligra",
  "//ligra:speculative_for",
  "//ligra/pbbslib:dyn_arr",
  "//ligra/pbbslib:sparse_table",
  "//pbbslib:random",
  "//pbbslib:sample_sort",
  "//pbbslib/strings:string_basics",
  ]
)

gbbs_cc_binary(
  name = "Biconnectivity_main",
  srcs = ["Biconnectivity.cc"],
  deps = [
  ":Biconnectivity",
  "//ligra:ligra",
  "//ligra/pbbslib:sparse_additive_map",
  "//pbbslib/strings:string_basics",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
