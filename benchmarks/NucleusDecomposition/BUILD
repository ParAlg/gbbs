cc_library(
  name = "commontable",
  hdrs = ["commontable.h"],
  srcs = ["commontable.cc"],
  deps = [
  "//benchmarks/CliqueCounting:relabel",
  "//benchmarks/CliqueCounting:induced_neighborhood",
  "//benchmarks/CliqueCounting:induced_intersection",
  "//benchmarks/CliqueCounting:induced_hybrid",
  "//benchmarks/CliqueCounting:induced_split",
  "//benchmarks/CliqueCounting:intersect",
  "//gbbs:bucket",
  "//gbbs:edge_map_reduce",
  "//gbbs:gbbs",
  "//gbbs/pbbslib:dyn_arr",
  "//gbbs/pbbslib:sparse_additive_map",
  "//gbbs/pbbslib:sparse_table",
  "//pbbslib:assert",
  "//pbbslib:list_allocator",
  "//pbbslib:integer_sort",
  ]
)

cc_library(
  name = "onetable",
  hdrs = ["onetable.h"],
  deps = [
  "//benchmarks/CliqueCounting:relabel",
  "//benchmarks/CliqueCounting:induced_neighborhood",
  "//benchmarks/CliqueCounting:induced_intersection",
  "//benchmarks/CliqueCounting:induced_hybrid",
  "//benchmarks/CliqueCounting:induced_split",
  "//benchmarks/CliqueCounting:intersect",
  "//gbbs:bucket",
  "//gbbs:edge_map_reduce",
  "//gbbs:gbbs",
  "//gbbs/pbbslib:dyn_arr",
  "//gbbs/pbbslib:sparse_additive_map",
  "//pbbslib:assert",
  "//pbbslib:list_allocator",
  "//pbbslib:integer_sort",
  ":commontable",
  ]
)

cc_library(
  name = "twotable",
  hdrs = ["twotable.h"],
  deps = [
  "//benchmarks/CliqueCounting:relabel",
  "//benchmarks/CliqueCounting:induced_neighborhood",
  "//benchmarks/CliqueCounting:induced_intersection",
  "//benchmarks/CliqueCounting:induced_hybrid",
  "//benchmarks/CliqueCounting:induced_split",
  "//benchmarks/CliqueCounting:intersect",
  "//gbbs:bucket",
  "//gbbs:edge_map_reduce",
  "//gbbs:gbbs",
  "//gbbs/pbbslib:dyn_arr",
  "//gbbs/pbbslib:sparse_additive_map",
  "//pbbslib:assert",
  "//pbbslib:list_allocator",
  "//pbbslib:integer_sort",
  ":commontable",
  ]
)

cc_library(
  name = "twotable_nosearch",
  hdrs = ["twotable_nosearch.h"],
  deps = [
  "//benchmarks/CliqueCounting:relabel",
  "//benchmarks/CliqueCounting:induced_neighborhood",
  "//benchmarks/CliqueCounting:induced_intersection",
  "//benchmarks/CliqueCounting:induced_hybrid",
  "//benchmarks/CliqueCounting:induced_split",
  "//benchmarks/CliqueCounting:intersect",
  "//gbbs:bucket",
  "//gbbs:edge_map_reduce",
  "//gbbs:gbbs",
  "//gbbs/pbbslib:dyn_arr",
  "//gbbs/pbbslib:sparse_additive_map",
  "//pbbslib:assert",
  "//pbbslib:list_allocator",
  "//pbbslib:integer_sort",
  ":commontable",
  ]
)

cc_library(
  name = "multitable",
  hdrs = ["multitable.h"],
  deps = [
  "//benchmarks/CliqueCounting:relabel",
  "//benchmarks/CliqueCounting:induced_neighborhood",
  "//benchmarks/CliqueCounting:induced_intersection",
  "//benchmarks/CliqueCounting:induced_hybrid",
  "//benchmarks/CliqueCounting:induced_split",
  "//benchmarks/CliqueCounting:intersect",
  "//gbbs:bucket",
  "//gbbs:edge_map_reduce",
  "//gbbs:gbbs",
  "//gbbs/pbbslib:dyn_arr",
  "//gbbs/pbbslib:sparse_additive_map",
  "//pbbslib:assert",
  "//pbbslib:list_allocator",
  "//pbbslib:integer_sort",
  ":commontable",
  ]
)

cc_library(
  name = "multitable_nosearch",
  hdrs = ["multitable_nosearch.h"],
  deps = [
  "//benchmarks/CliqueCounting:relabel",
  "//benchmarks/CliqueCounting:induced_neighborhood",
  "//benchmarks/CliqueCounting:induced_intersection",
  "//benchmarks/CliqueCounting:induced_hybrid",
  "//benchmarks/CliqueCounting:induced_split",
  "//benchmarks/CliqueCounting:intersect",
  "//gbbs:bucket",
  "//gbbs:edge_map_reduce",
  "//gbbs:gbbs",
  "//gbbs/pbbslib:dyn_arr",
  "//gbbs/pbbslib:sparse_additive_map",
  "//pbbslib:assert",
  "//pbbslib:list_allocator",
  "//pbbslib:integer_sort",
  ":commontable",
  ]
)

cc_library(
  name = "NucleusDecomposition",
  hdrs = ["NucleusDecomposition.h"],
  deps = [
  "//benchmarks/CliqueCounting:relabel",
  "//benchmarks/CliqueCounting:induced_neighborhood",
  "//benchmarks/CliqueCounting:induced_intersection",
  "//benchmarks/CliqueCounting:induced_hybrid",
  "//benchmarks/CliqueCounting:induced_split",
  "//benchmarks/CliqueCounting:intersect",
  "//gbbs:bucket",
  "//gbbs:edge_map_reduce",
  "//gbbs:gbbs",
  "//gbbs/pbbslib:dyn_arr",
  "//gbbs/pbbslib:sparse_additive_map",
  "//pbbslib:assert",
  "//pbbslib:list_allocator",
  "//pbbslib:integer_sort",
  "//benchmarks/DegeneracyOrder/BarenboimElkin08:DegeneracyOrder",
  "//benchmarks/DegeneracyOrder/GoodrichPszona11:DegeneracyOrder",
  "//benchmarks/KCore/JulienneDBS17:KCore",
  "//benchmarks/TriangleCounting/ShunTangwongsan15:Triangle",
  "//benchmarks/CliqueCounting:Clique",
  ":onetable",
  ":twotable",
  ":twotable_nosearch",
  ":multitable",
  ":multitable_nosearch",
  ":commontable",
  ]
)

cc_binary(
  name = "NucleusDecomposition_main",
  srcs = ["NucleusDecomposition.cc"],
  deps = [
  ":NucleusDecomposition"
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
