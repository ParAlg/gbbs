load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "GraphColoring",
  hdrs = ["GraphColoring.h"],
  deps = [
  "//ligra:ligra",
  "//pbbslib:random_shuffle"
  ]
)

gbbs_cc_binary(
  name = "GraphColoring_main",
  srcs = ["GraphColoring.cc"],
  deps = [":GraphColoring"]
)

package(
  default_visibility = ["//visibility:public"],
)
