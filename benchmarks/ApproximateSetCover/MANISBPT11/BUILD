cc_library(
  name = "ApproximateSetCover",
  hdrs = ["ApproximateSetCover.h"],
  deps = [
  "//gbbs:gbbs",
  "//gbbs:julienne",
  "//pbbslib:random",
  "//pbbslib:random_shuffle",
  ]
)

cc_binary(
  name = "ApproximateSetCover_main",
  srcs = ["ApproximateSetCover.cc"],
  deps = [":ApproximateSetCover"]
)
