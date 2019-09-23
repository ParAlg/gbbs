cc_library(
  name = "MinimumSpanningForest",
  hdrs = ["MinimumSpanningForest.h"],
  deps = [
  "//ligra:ligra",
  "//ligra:speculative_for",
  "//ligra:union_find",
  "//ligra/pbbslib:dyn_arr",
  "//pbbslib:binary_search",
  "//pbbslib:random",
  "//pbbslib:sample_sort",
  ]
)

cc_binary(
  name = "MinimumSpanningForest_main",
  srcs = ["MinimumSpanningForest.cpp"],
  deps = [":MinimumSpanningForest"]
)

package(
  default_visibility = ["//visibility:public"],
)
