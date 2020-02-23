cc_library(
  name = "intersect",
  hdrs = ["intersect.h"],
  deps = [
  "//ligra:bucket",
  "//ligra:edge_map_reduce",
  "//ligra:ligra",
  "//ligra/pbbslib:dyn_arr",
  "//ligra/pbbslib:sparse_table",
  "@simdinter//:intersection",
  "@graphsetinter//:intersection",
  ]
)

cc_library(
  name = "induced_intersection",
  hdrs = ["induced_intersection.h"],
  deps = [
  ":intersect",
  "//pbbslib:seq",
  ]
)

cc_library(
  name = "induced_neighborhood",
  hdrs = ["induced_neighborhood.h"],
  deps = [
  ":intersect",
  "//pbbslib:seq",
  ]
)

cc_library(
  name = "induced_hybrid",
  hdrs = ["induced_hybrid.h"],
  deps = [
  ":intersect",
  "//pbbslib:seq",
  ]
)

cc_library(
  name = "induced_split",
  hdrs = ["induced_split.h"],
  deps = [
  ":intersect",
  "//pbbslib:seq",
  ]
)

cc_library(
  name = "peel",
  hdrs = ["peel.h"],
  deps = [
  ":induced_hybrid",
  ":intersect",
  "//ligra:bucket",
  "//ligra:edge_map_reduce",
  "//ligra:ligra",
  "//ligra/pbbslib:dyn_arr",
  "//ligra/pbbslib:sparse_table",
  "//pbbslib:list_allocator",
  "//pbbslib:integer_sort",
  ]
)

cc_library(
  name = "relabel",
  hdrs = ["relabel.h"],
  deps = [
  "//ligra:ligra",
  "//ligra:graph",
  "//pbbslib:seq",
  ]
)

cc_library(
  name = "Clique",
  hdrs = ["Clique.h"],
  deps = [
  ":relabel",
  ":induced_neighborhood",
  ":induced_intersection",
  ":induced_hybrid",
  ":induced_split",
  ":peel",
  "//ligra:bucket",
  "//ligra:edge_map_reduce",
  "//ligra:ligra",
  "//ligra/pbbslib:dyn_arr",
  "//pbbslib:assert",
  "//pbbslib:list_allocator",
  "//pbbslib:integer_sort",
  "//benchmarks/DegeneracyOrder/BarenboimElkin08:DegeneracyOrder",
  "//benchmarks/DegeneracyOrder/GoodrichPszona11:DegeneracyOrder",
  "//benchmarks/KCore/JulienneDBS17:KCore",
  "//benchmarks/TriangleCounting/ShunTangwongsan15:Triangle",
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
