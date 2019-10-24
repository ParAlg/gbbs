cc_library(
  name = "intersect",
  hdrs = ["intersect.h"],
  deps = [
  "//ligra:bucket",
  "//ligra:edge_map_reduce",
  "//ligra:ligra",
  "//ligra/pbbslib:dyn_arr",
  "@simdinter//:intersection",
  "@graphsetinter//:intersection",
  ]
)

cc_library(
  name = "Clique",
  hdrs = ["Clique.h"],
  deps = [
  "//ligra:bucket",
  "//ligra:edge_map_reduce",
  "//ligra:ligra",
  "//ligra/pbbslib:dyn_arr",
  "//pbbslib:list_allocator",
  "//pbbslib:integer_sort",
  "//benchmarks/DegeneracyOrder/BarenboimElkin08:DegeneracyOrder",
  "//benchmarks/DegeneracyOrder/GoodrichPszona11:DegeneracyOrder",
  "//benchmarks/KCore/JulienneDBS17:KCore",
  ":intersect",
  ]
)

cc_binary(
  name = "Clique_main",
  srcs = ["Clique.cc"],
  deps = [
  ":Clique",
  "@simdinter//:intersectiondeps",],
  copts = ["-Iexternal/simdinter/include"],
)

package(
  default_visibility = ["//visibility:public"],
)
