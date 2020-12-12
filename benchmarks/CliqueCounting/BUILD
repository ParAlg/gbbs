cc_library(
  name = "intersect",
  hdrs = ["intersect.h"],
  deps = [
  "//gbbs:bucket",
  "//gbbs:edge_map_reduce",
  "//gbbs:gbbs",
  "//gbbs/pbbslib:dyn_arr",
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
  "//gbbs:bucket",
  "//gbbs:edge_map_reduce",
  "//gbbs:gbbs",
  "//gbbs/pbbslib:dyn_arr",
  "//gbbs/pbbslib:sparse_table",
  "//pbbslib:list_allocator",
  "//pbbslib:integer_sort",
  ]
)

cc_library(
  name = "relabel",
  hdrs = ["relabel.h"],
  deps = [
  "//gbbs:gbbs",
  "//gbbs:graph",
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
  "//gbbs:bucket",
  "//gbbs:edge_map_reduce",
  "//gbbs:gbbs",
  "//gbbs/pbbslib:dyn_arr",
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
  ":Clique"
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
