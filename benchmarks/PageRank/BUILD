cc_library(
  name = "PageRank",
  hdrs = ["PageRank.h"],
  deps = [
  "//gbbs:ligra",
  "//gbbs:edge_map_reduce",
  "//gbbs/pbbslib:sparse_table",
  ]
)

cc_binary(
  name = "PageRank_main",
  srcs = ["PageRank.cc"],
  deps = [":PageRank"]
)

package(
  default_visibility = ["//visibility:public"],
)
