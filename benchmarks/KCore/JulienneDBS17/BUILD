cc_library(
  name = "KCore",
  hdrs = ["KCore.h"],
  deps = [
  "//ligra:ligra",
  "//ligra:bucket"
  "//ligra:edge_map_reduce"
  "//ligra/pbbslib:dyn_arr"
  ]
)

cc_binary(
  name = "KCore_main",
  srcs = ["KCore.cpp"],
  deps = [":KCore"]
)

package(
  default_visibility = ["//visibility:public"],
)
