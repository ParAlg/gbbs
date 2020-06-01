cc_library(
  name = "MinimumSpanningForest",
  hdrs = ["MinimumSpanningForest.h"],
  deps = [
  "//gbbs:gbbs",
  "//gbbs:speculative_for",
  "//gbbs:union_find",
  "//gbbs/pbbslib:dyn_arr",
  "//pbbslib:binary_search",
  "//pbbslib:random",
  "//pbbslib:sample_sort",
  ]
)

cc_binary(
  name = "MinimumSpanningForest_main",
  srcs = ["MinimumSpanningForest.cc"],
  deps = [":MinimumSpanningForest"]
)

package(
  default_visibility = ["//visibility:public"],
)
