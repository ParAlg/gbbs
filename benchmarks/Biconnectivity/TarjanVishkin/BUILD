cc_library(
  name = "Biconnectivity",
  hdrs = ["Biconnectivity.h"],
  deps = [
  "//ligra:ligra"
  "//benchmarks/Connectivity/SDB13:Connectivity"
  "//ligra/pbbslib:dyn_arr",
  "//pbbslib:random",
  "//pbbslib:sample_sort",
  "//pbbslib:sparse_table",
  "//pbbslib:speculative_for",
  ]
)

cc_binary(
  name = "Biconnectivity_main",
  srcs = ["Biconnectivity.cpp"],
  deps = [":Biconnectivity"]
)
