cc_library(
  name = "Biconnectivity",
  hdrs = ["Biconnectivity.h"],
  deps = [
  "//benchmarks/Connectivity/WorkEfficientSDB14:Connectivity",
  "//ligra:gbbs",
  "//ligra/pbbslib:dyn_arr",
  "//ligra/pbbslib:sparse_table",
  "//pbbslib:random",
  "//pbbslib:sample_sort",
  "//pbbslib/strings:string_basics",
  ]
)

cc_binary(
  name = "Biconnectivity_main",
  srcs = ["Biconnectivity.cc"],
  deps = [
  ":Biconnectivity",
  "//ligra:ligra",
  "//ligra/pbbslib:sparse_additive_map",
  "//pbbslib/strings:string_basics",
  ]
)
