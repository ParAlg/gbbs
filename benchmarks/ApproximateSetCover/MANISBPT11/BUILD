cc_library(
  name = "ApproximateSetCover",
  hdrs = ["ApproximateSetCover.h"],
  deps = [
  "//ligra:gbbs",
  "//ligra:julienne",
  "//pbbslib:random",
  "//pbbslib:random_shuffle",
  ]
)

cc_binary(
  name = "ApproximateSetCover_main",
  srcs = ["ApproximateSetCover.cc"],
  deps = [":ApproximateSetCover"]
)
