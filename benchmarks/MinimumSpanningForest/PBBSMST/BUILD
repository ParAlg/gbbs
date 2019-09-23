cc_library(
  name = "MinimumSpanningForest",
  hdrs = ["MinimumSpanningForest.h"],
  deps = [
  "//ligra:ligra",
  "//ligra/pbbslib:sparse_table"
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
