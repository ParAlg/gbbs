cc_library(
  name = "Spanner",
  hdrs = ["Spanner.h"],
  deps = [
  "//ligra:ligra",
  "//ligra/pbbslib:sparse_table",
  "//ligra/pbbslib:dyn_arr",
  "//benchmarks/LowDiameterDecomposition/MPX13:LowDiameterDecomposition",
  ]
)

cc_binary(
  name = "Spanner_main",
  srcs = ["Spanner.cpp"],
  deps = [":Spanner"]
)

package(
  default_visibility = ["//visibility:public"],
)
