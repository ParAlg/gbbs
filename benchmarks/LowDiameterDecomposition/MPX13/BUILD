load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "LowDiameterDecomposition",
  hdrs = ["LowDiameterDecomposition.h"],
  deps = [
  "//pbbslib:random_shuffle",
  "//ligra:ligra"
  ]
)

gbbs_cc_binary(
  name = "LowDiameterDecomposition_main",
  srcs = ["LowDiameterDecomposition.cc"],
  deps = [":LowDiameterDecomposition"]
)

package(
  default_visibility = ["//visibility:public"],
)
