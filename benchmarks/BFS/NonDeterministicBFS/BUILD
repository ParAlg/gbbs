load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "BFS",
  hdrs = ["BFS.h"],
  deps = ["//ligra:ligra"]
)

gbbs_cc_binary(
  name = "BFS_main",
  srcs = ["BFS.cc"],
  deps = [":BFS"]
)

package(
  default_visibility = ["//visibility:public"],
)
