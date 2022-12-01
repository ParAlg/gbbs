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
  "//gbbs/helpers:dyn_arr",
  "//gbbs/helpers:sparse_additive_map",
  "//gbbs/helpers:sparse_table",
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
  "//gbbs/helpers:dyn_arr",
  "//gbbs/helpers:sparse_additive_map",
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
  "//gbbs/helpers:dyn_arr",
  "//gbbs/helpers:sparse_additive_map",
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
  "//gbbs/helpers:dyn_arr",
  "//gbbs/helpers:sparse_additive_map",
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
  "//gbbs/helpers:dyn_arr",
  "//gbbs/helpers:sparse_additive_map",
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
  "//gbbs/helpers:dyn_arr",
  "//gbbs/helpers:sparse_additive_map",
  ":commontable",
  ]
)

cc_library(
  name = "list_buffer",
  hdrs = ["list_buffer.h"],
  deps = [
  "//gbbs:bucket",
  "//gbbs:edge_map_reduce",
  "//gbbs:gbbs",
  "//gbbs/helpers:dyn_arr",
  "//gbbs/helpers:sparse_additive_map",
  ]
)

cc_library(
  name = "NucleusDecomposition_structs",
  hdrs = ["NucleusDecomposition_structs.h"],
  deps = [
  "//benchmarks/CliqueCounting:relabel",
  "//benchmarks/CliqueCounting:induced_neighborhood",
  "//benchmarks/CliqueCounting:induced_intersection",
  "//benchmarks/CliqueCounting:induced_hybrid",
  "//benchmarks/CliqueCounting:induced_split",
  "//benchmarks/CliqueCounting:intersect",
  "//benchmarks/Connectivity/SimpleUnionAsync:Connectivity",
  "//gbbs:bucket",
  "//gbbs:edge_map_reduce",
  "//gbbs:gbbs",
  "//gbbs/helpers:dyn_arr",
  "//gbbs/helpers:sparse_additive_map",
  "//benchmarks/DegeneracyOrder/BarenboimElkin08:DegeneracyOrder",
  "//benchmarks/DegeneracyOrder/GoodrichPszona11:DegeneracyOrder",
  "//benchmarks/KCore/JulienneDBS17:KCore",
  "//benchmarks/TriangleCounting/ShunTangwongsan15:Triangle",
  "//benchmarks/CliqueCounting:Clique",
  ":twotable_nosearch",
  ":commontable",
  ":list_buffer",
  ]
)

cc_library(
  name = "NucleusDecomposition_pnd",
  hdrs = ["NucleusDecomposition_pnd.h"],
  deps = [
  "//benchmarks/CliqueCounting:relabel",
  "//benchmarks/CliqueCounting:induced_neighborhood",
  "//benchmarks/CliqueCounting:induced_intersection",
  "//benchmarks/CliqueCounting:induced_hybrid",
  "//benchmarks/CliqueCounting:induced_split",
  "//benchmarks/CliqueCounting:intersect",
  "//benchmarks/Connectivity/SimpleUnionAsync:Connectivity",
  "//gbbs:bucket",
  "//gbbs:edge_map_reduce",
  "//gbbs:gbbs",
  "//gbbs/helpers:dyn_arr",
  "//gbbs/helpers:sparse_additive_map",
  "//benchmarks/DegeneracyOrder/BarenboimElkin08:DegeneracyOrder",
  "//benchmarks/DegeneracyOrder/GoodrichPszona11:DegeneracyOrder",
  "//benchmarks/KCore/JulienneDBS17:KCore",
  "//benchmarks/TriangleCounting/ShunTangwongsan15:Triangle",
  "//benchmarks/CliqueCounting:Clique",
  ":twotable_nosearch",
  ":commontable",
  ":list_buffer",
  ]
)

cc_library(
  name = "NucleusDecomposition_common",
  hdrs = ["NucleusDecomposition_common.h"],
  deps = [
  "//benchmarks/CliqueCounting:relabel",
  "//benchmarks/CliqueCounting:induced_neighborhood",
  "//benchmarks/CliqueCounting:induced_intersection",
  "//benchmarks/CliqueCounting:induced_hybrid",
  "//benchmarks/CliqueCounting:induced_split",
  "//benchmarks/CliqueCounting:intersect",
  "//benchmarks/Connectivity/SimpleUnionAsync:Connectivity",
  "//gbbs:bucket",
  "//gbbs:edge_map_reduce",
  "//gbbs:gbbs",
  "//gbbs/helpers:dyn_arr",
  "//gbbs/helpers:sparse_additive_map",
  "//benchmarks/DegeneracyOrder/BarenboimElkin08:DegeneracyOrder",
  "//benchmarks/DegeneracyOrder/GoodrichPszona11:DegeneracyOrder",
  "//benchmarks/KCore/JulienneDBS17:KCore",
  "//benchmarks/TriangleCounting/ShunTangwongsan15:Triangle",
  "//benchmarks/CliqueCounting:Clique",
  ":twotable_nosearch",
  ":commontable",
  ":list_buffer",
  ":NucleusDecomposition_structs",
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
  "//gbbs",
  "//gbbs/helpers:dyn_arr",
  "//gbbs/helpers:sparse_additive_map",
  "//benchmarks/DegeneracyOrder/BarenboimElkin08:DegeneracyOrder",
  "//benchmarks/DegeneracyOrder/GoodrichPszona11:DegeneracyOrder",
  "//benchmarks/KCore/JulienneDBS17:KCore",
  "//benchmarks/TriangleCounting/ShunTangwongsan15:Triangle",
  "//benchmarks/CliqueCounting:Clique",
  ":twotable_nosearch",
  ":commontable",
  ":NucleusDecomposition_common",
  ]
)

cc_library(
  name = "TwoThreeNucleusDecomposition",
  hdrs = ["TwoThreeNucleusDecomposition.h"],
  deps = [
  "//benchmarks/CliqueCounting:relabel",
  "//benchmarks/CliqueCounting:induced_neighborhood",
  "//benchmarks/CliqueCounting:induced_intersection",
  "//benchmarks/CliqueCounting:induced_hybrid",
  "//benchmarks/CliqueCounting:induced_split",
  "//benchmarks/CliqueCounting:intersect",
  "//gbbs:bucket",
  "//gbbs:edge_map_reduce",
  "//gbbs",
  "//gbbs/helpers:dyn_arr",
  "//gbbs/helpers:sparse_additive_map",
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
  ":NucleusDecomposition_common",
  ]
)

cc_binary(
  name = "NucleusDecomposition_main",
  srcs = ["NucleusDecomposition.cc"],
  deps = [
  ":NucleusDecomposition",
  "//benchmarks/KTruss:KTruss",
  ]
)

cc_binary(
  name = "TwoThreeNucleusDecomposition_main",
  srcs = ["TwoThreeNucleusDecomposition.cc"],
  deps = [
  ":TwoThreeNucleusDecomposition"
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
