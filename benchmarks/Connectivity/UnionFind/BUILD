cc_library(
  name = "union_find_rules",
  hdrs = ["union_find_rules.h"],
  deps = [
  "//ligra:bridge",
  "//pbbslib:seq"
  ]
)

cc_library(
  name = "Connectivity",
  hdrs = ["Connectivity.h"],
  deps = [
  "//ligra:bridge",
  "//ligra:ligra",
  "//ligra/pbbslib:sparse_table"
  ]
)

cc_binary(
  name = "Connectivity_main",
  srcs = ["Connectivity.cpp"],
  deps = [
  ":Connectivity",
  ":union_find_rules",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
