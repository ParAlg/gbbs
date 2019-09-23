cc_library(
  name = "BellmanFord",
  hdrs = ["BellmanFord.h"],
  deps = [
  "//ligra:ligra",
  ]
)

cc_binary(
  name = "BellmanFord_main",
  srcs = ["BellmanFord.cpp"],
  deps = [":BellmanFord"]
)

package(
  default_visibility = ["//visibility:public"],
)
