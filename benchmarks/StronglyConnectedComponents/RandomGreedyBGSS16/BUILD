load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "StronglyConnectedComponents",
  hdrs = ["StronglyConnectedComponents.h"],
  deps = [
  "//ligra:ligra",
  "//ligra/pbbslib:resizable_table",
  "//pbbslib:random_shuffle"
  ]
)

gbbs_cc_binary(
  name = "StronglyConnectedComponents_main",
  srcs = ["StronglyConnectedComponents.cc"],
  deps = [":StronglyConnectedComponents"]
)

package(
  default_visibility = ["//visibility:public"],
)
