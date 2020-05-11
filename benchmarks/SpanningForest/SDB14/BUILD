cc_library(
  name = "contract_sf",
  hdrs = ["contract_sf.h"],
  deps = [
  "//gbbs:graph",
  "//gbbs:macros",
  "//gbbs/pbbslib:sparse_table",
  "//gbbs/pbbslib:dyn_arr",
  "//pbbslib:sequence_ops"
  ]
)

cc_library(
  name = "SpanningForest",
  hdrs = ["SpanningForest.h"],
  deps = [
  ":contract_sf",
  "//benchmarks/LowDiameterDecomposition/MPX13:LowDiameterDecomposition",
  "//gbbs:ligra",
  "//gbbs/pbbslib:sparse_table"
  ]
)

cc_binary(
  name = "SpanningForest_main",
  srcs = ["SpanningForest.cc"],
  deps = [
  ":SpanningForest",
  "//benchmarks/SpanningForest:common",
  "//benchmarks/SpanningForest/BFSSF:SpanningForest"
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
