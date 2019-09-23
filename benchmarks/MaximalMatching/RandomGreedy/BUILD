cc_library(
  name = "MaximalMatching",
  hdrs = ["MaximalMatching.h"],
  deps = [
  "//ligra:bridge",
  "//ligra:ligra",
  "//ligra:speculative_for",
  "//ligra/pbbslib:sparse_table",
  "//ligra/pbbslib:dyn_arr"
  ]
)

cc_binary(
  name = "MaximalMatching_main",
  srcs = ["MaximalMatching.cpp"],
  deps = [
  ":MaximalMatching",
  "//ligra:bridge",
  "//ligra:ligra",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
