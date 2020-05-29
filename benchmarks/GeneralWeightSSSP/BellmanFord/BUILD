load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "BellmanFord",
  hdrs = ["BellmanFord.h"],
  deps = [
  "//ligra:ligra",
  ]
)

gbbs_cc_binary(
  name = "BellmanFord_main",
  srcs = ["BellmanFord.cc"],
  deps = [":BellmanFord"]
)

package(
  default_visibility = ["//visibility:public"],
)
