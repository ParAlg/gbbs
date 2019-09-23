cc_library(
  name = "SSBetweennessCentrality",
  hdrs = ["SSBetweennessCentrality.h"],
  deps = [
  "//ligra:ligra",
  "//ligra/pbbslib:sparse_table"
  "//ligra/pbbslib:edge_map_reduce"
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
