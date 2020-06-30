cc_library(
  name = "Connectivity",
  hdrs = ["Connectivity.h"],
  deps = [
  "//benchmarks/Connectivity:common",
  "//gbbs:gbbs",
  ]
)

cc_binary(
  name = "Connectivity_main",
  srcs = ["Connectivity.cc"],
  deps = [":Connectivity"]
)

package(
  default_visibility = ["//visibility:public"],
)
