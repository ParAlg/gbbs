cc_library(
  name = "MaximalMatching",
  hdrs = ["MaximalMatching.h"],
  deps = [
  "//gbbs:gbbs",
  "//gbbs/helpers:sparse_table",
  "//gbbs/helpers:speculative_for",
  "//gbbs/helpers:dyn_arr"
  ]
)

cc_binary(
  name = "MaximalMatching_main",
  srcs = ["MaximalMatching.cc"],
  deps = [
  ":MaximalMatching",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
