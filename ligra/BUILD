load("//internal_tools:build_defs.bzl", "gbbs_cc_library")

gbbs_cc_library(
  name = "bridge",
  hdrs = ["bridge.h"],
  srcs = ["bridge.cc"],
  deps = [
  "//pbbslib:binary_search",
  "//pbbslib:counting_sort",
  "//pbbslib:integer_sort",
  "//pbbslib:monoid",
  "//pbbslib:parallel",
  "//pbbslib:random",
  "//pbbslib:random_shuffle",
  "//pbbslib:sample_sort",
  "//pbbslib:seq",
  "//pbbslib:sequence_ops",
  "//pbbslib:utilities",
  ]
)

gbbs_cc_library(
  name = "io",
  hdrs = ["io.h"],
  srcs = ["io.cc"],
  deps = [
  ":macros",
  ]
)

gbbs_cc_library(
  name = "graph_io",
  hdrs = ["graph_io.h"],
  srcs = ["graph_io.cc"],
  deps = [
  ":graph",
  ":io",
  ":macros",
  ":vertex",
  "//pbbslib:sample_sort",
  "//pbbslib:stlalgs",
  "//pbbslib/strings:string_basics",
  ]
)

gbbs_cc_library(
  name = "bucket",
  hdrs = ["bucket.h"],
  deps = [
  ":vertex_subset",
  ":bridge",
  "//ligra/pbbslib:dyn_arr"
  ]
)

gbbs_cc_library(
  name = "compressed_vertex",
  hdrs = ["compressed_vertex.h"],
  deps = [
  ":macros",
  "//pbbslib:monoid",
  "//ligra/encodings:decoders"
  ]
)

gbbs_cc_library(
  name = "edge_map_blocked",
  hdrs = ["edge_map_blocked.h"],
  srcs = ["edge_map_blocked.cc"],
  deps = [
  ":bridge",
  ":edge_map_utils",
  ":graph",
  ":vertex_subset",
  "//pbbslib:binary_search",
  "//pbbslib:list_allocator",
  ]
)

gbbs_cc_library(
  name = "edge_map_reduce",
  hdrs = ["edge_map_reduce.h"],
  deps = [
  ":bridge",
  ":histogram",
  ":ligra",
  ":flags",
  ]
)

gbbs_cc_library(
  name = "edge_map_utils",
  hdrs = ["edge_map_utils.h"],
  deps = [
  ":macros",
  "//pbbslib:binary_search",
  "//pbbslib:utilities",
  ]
)

gbbs_cc_library(
  name = "flags",
  hdrs = ["flags.h"],
)

gbbs_cc_library(
  name = "graph",
  hdrs = ["graph.h"],
  srcs = ["graph.cc"],
  deps = [
  ":bridge",
  ":compressed_vertex",
  ":flags",
  ":vertex",
  ]
)

gbbs_cc_library(
  name = "graph_mutation",
  hdrs = ["graph_mutation.h"],
  deps = [
  ":bridge",
  ":graph",
  ]
)

gbbs_cc_library(
  name = "graph_test_utils",
  hdrs = ["graph_test_utils.h"],
  srcs = ["graph_test_utils.cc"],
  deps = [
  ":graph",
  ":macros",
  ":undirected_edge",
  ":vertex",
  "//pbbslib:utilities",
  "@googletest//:gtest_main",
  ]
)

gbbs_cc_library(
  name = "histogram",
  hdrs = ["histogram.h"],
  deps = [
  ":bridge",
  ":macros",
  ":sequential_ht",
  "//ligra/pbbslib:counting_sort_no_transpose",
  ]
)

gbbs_cc_library(
  name = "ligra",
  hdrs = ["ligra.h"],
  srcs = ["ligra.cc"],
  deps = [
  ":bridge",
  ":compressed_vertex",
  ":edge_map_utils",
  ":edge_map_blocked",
  ":flags",
  ":graph",
  ":graph_io",
  ":graph_mutation",
  ":parse_command_line",
  ":vertex",
  ":vertex_subset",
  ]
)

gbbs_cc_library(
  name = "macros",
  hdrs = ["macros.h"],
  deps = [
  ":bridge",
  ]
)

gbbs_cc_library(
  name = "parse_command_line",
  hdrs = ["parse_command_line.h"],
  srcs = ["parse_command_line.cc"],
)

gbbs_cc_library(
  name = "sequential_ht",
  hdrs = ["sequential_ht.h"],
  deps = [
  ":macros",
  "//pbbslib:utilities"
  ]
)

gbbs_cc_library(
  name = "speculative_for",
  hdrs = ["speculative_for.h"],
  deps = [
  ":bridge"
  ]
)

gbbs_cc_library(
  name = "union_find",
  hdrs = ["union_find.h"],
  srcs = ["union_find.cc"],
  deps = [
  ":macros",
  ":speculative_for"
  ]
)

gbbs_cc_library(
  name = "vertex",
  hdrs = ["vertex.h"],
  deps = [
  ":macros",
  ":speculative_for"
  ]
)

gbbs_cc_library(
  name = "vertex_subset",
  hdrs = ["vertex_subset.h"],
  deps = [
  ":bridge",
  ":macros",
  ]
)

gbbs_cc_library(
  name = "undirected_edge",
  hdrs = ["undirected_edge.h"],
  srcs = ["undirected_edge.cc"],
  deps = [
  ":macros",
  "//pbbslib:utilities"
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
