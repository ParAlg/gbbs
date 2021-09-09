licenses(["notice"])

cc_library(
  name = "Spanner",
  hdrs = ["Spanner.h"],
  deps = [
  "//gbbs:gbbs",
  "//gbbs/helpers:sparse_table",
  "//gbbs/helpers:dyn_arr",
  "//benchmarks/LowDiameterDecomposition/MPX13:LowDiameterDecomposition",
  ]
)

cc_binary(
  name = "Spanner_main",
  srcs = ["Spanner.cc"],
  deps = [":Spanner"]
)

package(
  default_visibility = ["//visibility:public"],
)
