cc_library(
  name = "Biconnectivity",
  hdrs = ["Biconnectivity.h"],
  deps = [
  "//benchmarks/Connectivity/WorkEfficientSDB14:Connectivity",
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
  srcs = ["Biconnectivity.cc"],
  deps = [
  ":Biconnectivity",
  "//ligra/pbbslib:sparse_additive_map",
  "//pbbslib/strings:string_basics",
  ]
)
