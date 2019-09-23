cc_library(
  name = "SSBetweennessCentrality",
  hdrs = ["SSBetweennessCentrality.h"],
  deps = [
  "//ligra:ligra",
  "//ligra:edge_map_reduce",
  "//ligra/pbbslib:sparse_table",
  ]
)

cc_binary(
  name = "SSBetweennessCentrality_main",
  srcs = ["SSBetweennessCentrality.cpp"],
  deps = [":SSBetweennessCentrality"]
)

package(
  default_visibility = ["//visibility:public"],
)
