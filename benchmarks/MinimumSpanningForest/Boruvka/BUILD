load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "MinimumSpanningForest",
  hdrs = ["MinimumSpanningForest.h"],
  deps = [
  "//ligra:ligra",
  "//ligra:union_find",
  "//ligra/pbbslib:sparse_table",
  "//ligra/pbbslib:dyn_arr",
  "//pbbslib:binary_search",
  "//pbbslib:random",
  "//pbbslib:sample_sort",
  ]
)

gbbs_cc_binary(
  name = "MinimumSpanningForest_main",
  srcs = ["MinimumSpanningForest.cc"],
  deps = [":MinimumSpanningForest"]
)

package(
  default_visibility = ["//visibility:public"],
)
