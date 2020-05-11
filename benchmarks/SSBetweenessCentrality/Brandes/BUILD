cc_library(
  name = "SSBetweennessCentrality",
  hdrs = ["SSBetweennessCentrality.h"],
  deps = [
  "//gbbs:gbbs",
  "//gbbs:edge_map_reduce",
  "//gbbs/pbbslib:sparse_table",
  ]
)

cc_binary(
  name = "SSBetweennessCentrality_main",
  srcs = ["SSBetweennessCentrality.cc"],
  deps = [":SSBetweennessCentrality"]
)

package(
  default_visibility = ["//visibility:public"],
)
