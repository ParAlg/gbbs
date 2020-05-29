load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "CoSimRank",
  hdrs = ["CoSimRank.h"],
  deps = [
  "//ligra:ligra",
  "//ligra:edge_map_reduce",
  ]
)

gbbs_cc_binary(
  name = "CoSimRank_main",
  srcs = ["CoSimRank.cc"],
  deps = [":CoSimRank"]
)

package(
  default_visibility = ["//visibility:public"],
)
