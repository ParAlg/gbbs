cc_library(
  name = "MinimumSpanningForest",
  hdrs = ["MinimumSpanningForest.h"],
  deps = [
  "//gbbs:gbbs",
  "//gbbs:union_find",
  "//gbbs/pbbslib:sparse_table",
  "//gbbs/pbbslib:dyn_arr",
  "//benchmarks/SpanningForest/SDB14:SpanningForest",
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
