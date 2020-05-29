load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "DeltaStepping",
  hdrs = ["DeltaStepping.h"],
  deps = [
  "//ligra:ligra",
  "//ligra:bucket",
  ]
)

gbbs_cc_binary(
  name = "DeltaStepping_main",
  srcs = ["DeltaStepping.cc"],
  deps = [":DeltaStepping"]
)

package(
  default_visibility = ["//visibility:public"],
)
