cc_library(
  name = "Biconnectivity",
  hdrs = ["Biconnectivity.h"],
  deps = [
  "//benchmarks/Connectivity/WorkEfficientSDB14:Connectivity",
  "//gbbs:gbbs",
  "//gbbs/pbbslib:dyn_arr",
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
  "//gbbs/pbbslib:sparse_additive_map",
  "//pbbslib/strings:string_basics",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
