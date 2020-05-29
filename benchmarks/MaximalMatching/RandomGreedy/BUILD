load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
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

gbbs_cc_binary(
  name = "MaximalMatching_main",
  srcs = ["MaximalMatching.cc"],
  deps = [
  ":MaximalMatching",
  "//ligra:bridge",
  "//ligra:ligra",
  "//pbbslib/strings:string_basics"
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
