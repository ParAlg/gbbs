cc_library(
  name = "MaximalMatching",
  hdrs = ["MaximalMatching.h"],
  deps = [
  "//ligra:gbbs",
  "//ligra:speculative_for",
  "//ligra/pbbslib:sparse_table",
  "//ligra/pbbslib:dyn_arr"
  ]
)

cc_binary(
  name = "MaximalMatching_main",
  srcs = ["MaximalMatching.cc"],
  deps = [
  ":MaximalMatching",
  "//ligra:gbbs",
  "//pbbslib/strings:string_basics"
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
