cc_library(
  name = "contract_sf",
  hdrs = ["contract_sf.h"],
  deps = [
  "//ligra:graph",
  "//ligra:macros",
  "//ligra/pbbslib:sparse_table",
  "//pbbslib:sequence_ops"
  ]
)

cc_library(
  name = "SpanningForest",
  hdrs = ["SpanningForest.h"],
  deps = [
  ":contract",
  "//benchmarks/LowDiameterDecomposition/MPX13:LowDiameterDecomposition",
  "//ligra:ligra",
  "//ligra/pbbslib:sparse_table"
  ]
)

cc_binary(
  name = "SpanningForest_main",
  srcs = ["SpanningForest.cpp"],
  deps = [":SpanningForest"]
)

package(
  default_visibility = ["//visibility:public"],
)
