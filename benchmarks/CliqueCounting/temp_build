cc_library(
  name = "intersect",
  hdrs = ["intersect.h"],
  deps = [
  "//ligra:bucket",
  "//ligra:edge_map_reduce",
  "//ligra:ligra",
  "//ligra/pbbslib:dyn_arr",
  "@simdinter//:intersection"
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
  "//benchmarks/ApproximateDensestSubgraph/GreedyCharikar:DensestSubgraph",
  "//benchmarks/ApproximateDensestSubgraph/ApproxPeelingBKV12:DensestSubgraph",
  ":intersect"
  ]
)

cc_binary(
  name = "Clique_main",
  srcs = ["Clique.cc"],
  deps = [":Clique"]
)

package(
  default_visibility = ["//visibility:public"],
)
