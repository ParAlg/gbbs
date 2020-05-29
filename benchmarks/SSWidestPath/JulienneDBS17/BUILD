load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "SSWidestPath",
  hdrs = ["SSWidestPath.h"],
  deps = [
  "//ligra:ligra",
  "//ligra:bucket"
  ]
)

gbbs_cc_binary(
  name = "SSWidestPath_main",
  srcs = ["SSWidestPath.cc"],
  deps = [":SSWidestPath"]
)

package(
  default_visibility = ["//visibility:public"],
)
