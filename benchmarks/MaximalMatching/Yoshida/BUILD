load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "reorder",
  hdrs = ["reorder.h"],
  deps = [
  "//ligra:bridge",
  "//ligra:ligra",
  ]
)

gbbs_cc_library(
  name = "yoshida_matching_utils",
  hdrs = ["yoshida_matching_utils.h"],
  deps = []
)

gbbs_cc_library(
  name = "MaximalMatching",
  hdrs = ["MaximalMatching.h"],
  deps = [
  ":reorder",
  ":yoshida_matching_utils",
  "//ligra:bridge",
  "//ligra:ligra",
  ]
)

gbbs_cc_binary(
  name = "MaximalMatching_main",
  srcs = ["MaximalMatching.cc"],
  deps = [
  ":MaximalMatching",
  "//ligra:bridge",
  "//ligra:ligra",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
