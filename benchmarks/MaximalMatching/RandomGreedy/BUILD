cc_library(
  name = "MaximalMatching",
  hdrs = ["MaximalMatching.h"],
  deps = [
  "//gbbs:gbbs",
  "//gbbs:speculative_for",
  "//gbbs/pbbslib:sparse_table",
  "//gbbs/pbbslib:dyn_arr"
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
