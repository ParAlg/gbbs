load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "KCore",
  hdrs = ["KCore.h"],
  deps = [
  "//ligra:ligra",
  "//ligra:bucket",
  "//ligra:edge_map_reduce",
  "//ligra/pbbslib:dyn_arr",
  ]
)

gbbs_cc_binary(
  name = "KCore_main",
  srcs = ["KCore.cc"],
  deps = [":KCore"]
)

package(
  default_visibility = ["//visibility:public"],
)
