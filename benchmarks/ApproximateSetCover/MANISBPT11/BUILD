cc_library(
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

cc_binary(
  name = "ApproximateSetCover_main",
  srcs = ["ApproximateSetCover.cc"],
  deps = [":ApproximateSetCover"]
)
