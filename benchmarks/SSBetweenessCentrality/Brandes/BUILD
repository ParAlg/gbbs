load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "SSBetweennessCentrality",
  hdrs = ["SSBetweennessCentrality.h"],
  deps = [
  "//ligra:ligra",
  "//ligra:edge_map_reduce",
  "//ligra/pbbslib:sparse_table",
  ]
)

gbbs_cc_binary(
  name = "SSBetweennessCentrality_main",
  srcs = ["SSBetweennessCentrality.cc"],
  deps = [":SSBetweennessCentrality"]
)

package(
  default_visibility = ["//visibility:public"],
)
