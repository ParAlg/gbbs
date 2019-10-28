cc_library(
  name = "contract",
  hdrs = ["contract.h"],
  deps = [
  "//ligra:graph",
  "//ligra:macros",
  "//ligra/pbbslib:sparse_table",
  "//pbbslib:sequence_ops"
  ]
)

cc_library(
  name = "Connectivity",
  hdrs = ["Connectivity.h"],
  deps = [
  ":contract",
  "//benchmarks/LowDiameterDecomposition/MPX13:LowDiameterDecomposition",
  "//ligra:ligra",
  "//ligra/pbbslib:sparse_table"
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
