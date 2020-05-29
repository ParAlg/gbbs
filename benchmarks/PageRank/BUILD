load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "PageRank",
  hdrs = ["PageRank.h"],
  deps = [
  "//ligra:ligra",
  "//ligra:edge_map_reduce",
  "//ligra/pbbslib:sparse_table",
  ]
)

gbbs_cc_binary(
  name = "PageRank_main",
  srcs = ["PageRank.cc"],
  deps = [":PageRank"]
)

package(
  default_visibility = ["//visibility:public"],
)
