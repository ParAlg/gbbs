load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "wBFS",
  hdrs = ["wBFS.h"],
  deps = [
  "//ligra:ligra",
  "//ligra:bucket",
  ]
)

gbbs_cc_binary(
  name = "wBFS_main",
  srcs = ["wBFS.cc"],
  deps = [":wBFS"]
)

package(
  default_visibility = ["//visibility:public"],
)
