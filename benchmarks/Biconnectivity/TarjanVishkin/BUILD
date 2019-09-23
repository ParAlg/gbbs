cc_library(
  name = "Biconnectivity",
  hdrs = ["Biconnectivity.h"],
  deps = [
  "//benchmarks/Connectivity/SDB14:Connectivity",
  "//ligra:ligra",
  "//ligra:speculative_for",
  "//ligra/pbbslib:dyn_arr",
  "//ligra/pbbslib:sparse_table",
  "//pbbslib:random",
  "//pbbslib:sample_sort",
  ]
)

cc_binary(
  name = "Biconnectivity_main",
  srcs = ["Biconnectivity.cpp"],
  deps = [
  ":Biconnectivity",
  "//ligra/pbbslib:sparse_additive_map",
  ]
)
