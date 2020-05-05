package(
  default_visibility = ["//visibility:public"],
)

cc_library(
  name = "graph_filter",
  hdrs = ["graph_filter.h"],
  deps = [
  ":bitset_managers",
  ":block_vertex",
  "//ligra:graph",
  ]
)

cc_library(
  name = "bitset_managers",
  hdrs = ["bitset_managers.h"],
  deps = [
  ":bitset",
  "//ligra:flags",
  "//ligra:macros",
  "//ligra/encodings:byte_pd_amortized",
  ]
)

cc_library(
  name = "block_vertex",
  hdrs = ["block_vertex.h"],
  deps = [
  ":bitset_managers",
  "//ligra:macros",
  "//pbbslib:sequence_ops",
  ]
)

cc_library(
  name = "bitset",
  hdrs = ["bitset.h"],
  deps = [
    "//ligra:macros"
  ]
)
