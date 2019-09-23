cc_library(
  name = "SSWidestPath",
  hdrs = ["SSWidestPath.h"],
  deps = [
  "//ligra:ligra",
  "//ligra:bucket"
  ]
)

cc_binary(
  name = "SSWidestPath_main",
  srcs = ["SSWidestPath.cpp"],
  deps = [":SSWidestPath"]
)

package(
  default_visibility = ["//visibility:public"],
)
