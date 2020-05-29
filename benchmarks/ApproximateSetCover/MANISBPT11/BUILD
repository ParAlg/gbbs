load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "ApproximateSetCover",
  hdrs = ["ApproximateSetCover.h"],
  deps = [
  "//ligra:ligra",
  "//ligra:bucket",
  "//ligra:edge_map_reduce",
  "//pbbslib:random",
  "//pbbslib:random_shuffle",
  ]
)

gbbs_cc_binary(
  name = "ApproximateSetCover_main",
  srcs = ["ApproximateSetCover.cc"],
  deps = [":ApproximateSetCover"]
)
