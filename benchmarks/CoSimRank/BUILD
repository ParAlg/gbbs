cc_library(
  name = "CoSimRank",
  hdrs = ["CoSimRank.h"],
  deps = [
  "//ligra:ligra",
  "//ligra:edge_map_reduce",
  ]
)

cc_binary(
  name = "CoSimRank_main",
  srcs = ["CoSimRank.cc"],
  deps = [":CoSimRank"]
)

package(
  default_visibility = ["//visibility:public"],
)
