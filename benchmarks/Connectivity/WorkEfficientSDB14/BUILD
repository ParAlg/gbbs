cc_library(
  name = "Connectivity",
  hdrs = ["Connectivity.h"],
  deps = [
  "//benchmarks/LowDiameterDecomposition/MPX13:LowDiameterDecomposition",
  "//benchmarks/Connectivity:common",
  "//gbbs:gbbs",
  "//gbbs:contract",
  "//gbbs/pbbslib:sparse_table"
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
