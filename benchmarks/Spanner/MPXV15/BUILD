cc_library(
  name = "Spanner",
  hdrs = ["Spanner.h"],
  deps = [
  "//gbbs:ligra",
  "//gbbs/pbbslib:sparse_table",
  "//gbbs/pbbslib:dyn_arr",
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
