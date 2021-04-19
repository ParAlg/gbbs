cc_library(
  name = "MinimumSpanningForest",
  hdrs = ["MinimumSpanningForest.h"],
  deps = [
        "//benchmarks/Connectivity/UnionFind:union_find_rules",
  "//gbbs:gbbs",
  "//gbbs:union_find",
  "//gbbs/pbbslib:sparse_table",
  "//gbbs/pbbslib:dyn_arr",
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
