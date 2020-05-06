cc_library(
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

cc_library(
  name = "io",
  hdrs = ["io.h"],
  srcs = ["io.cc"],
  deps = [
  ":macros",
  ]
)

cc_library(
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

cc_library(
  name = "bucket",
  hdrs = ["bucket.h"],
  deps = [
  ":vertex_subset",
  ":bridge",
  "//ligra/pbbslib:dyn_arr"
  ]
)

cc_library(
  name = "compressed_vertex",
  hdrs = ["compressed_vertex.h"],
  deps = [
  ":macros",
  "//pbbslib:monoid",
  "//ligra/encodings:decoders"
  ]
)

cc_library(
  name = "edge_map_blocked",
  hdrs = ["edge_map_blocked.h"],
  srcs = ["edge_map_blocked.cc"],
  deps = [
  ":bridge",
  ":edge_map_utils",
  ":flags",
  ":vertex_subset",
  "//pbbslib:binary_search",
  "//pbbslib:list_allocator",
  ]
)

cc_library(
  name = "edge_map_reduce",
  hdrs = ["edge_map_reduce.h"],
  deps = [
  ":bridge",
  ":histogram",
  ":ligra",
  ":flags",
  ]
)

cc_library(
  name = "edge_map_utils",
  hdrs = ["edge_map_utils.h"],
  deps = [
  ":macros",
  "//pbbslib:binary_search",
  "//pbbslib:utilities",
  ]
)

cc_library(
  name = "flags",
  hdrs = ["flags.h"],
)

cc_library(
  name = "edge_array",
  hdrs = ["edge_array.h"],
  deps = [
  ":bridge",
  ":macros",
  ":vertex",
  ":vertex_subset",
  ]
)

cc_library(
  name = "graph",
  hdrs = ["graph.h"],
  srcs = ["graph.cc"],
  deps = [
  ":bridge",
  ":compressed_vertex",
  ":edge_array",
  ":edge_map_reduce",
  ":flags",
  ":ligra",
  ":macros",
  ":vertex",
  ":vertex_subset",
  ]
)

cc_library(
  name = "graph_mutation",
  hdrs = ["graph_mutation.h"],
  deps = [
  ":bridge",
  ":graph",
  ]
)

cc_library(
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

cc_library(
  name = "histogram",
  hdrs = ["histogram.h"],
  deps = [
  ":bridge",
  ":macros",
  ":sequential_ht",
  "//ligra/pbbslib:counting_sort_no_transpose",
  ]
)

cc_library(
  name = "ligra",
  hdrs = ["ligra.h"],
  srcs = ["ligra.cc"],
  deps = [
  ":bridge",
  ":compressed_vertex",
  ":edge_map_utils",
  ":edge_map_blocked",
  ":flags",
  ":vertex",
  ":vertex_subset",
  ]
)

# main header
cc_library(
  name = "gbbs",
  hdrs = ["gbbs.h"],
  deps = [
  ":benchmark",
  ":bridge",
  ":graph",
  ":macros",
  ":vertex_subset",
  ]
)



# main header
cc_library(
  name = "julienne",
  hdrs = ["julienne.h"],
  deps = [
  ":gbbs",
  ":bucket",
  "//ligra/pbbslib:dyn_arr",
  ]
)



cc_library(
  name = "benchmark",
  hdrs = ["benchmark.h"],
  srcs = ["benchmark.cc"],
  deps = [
  ":graph_io",
  ":parse_command_line",
  ]
)

cc_library(
  name = "macros",
  hdrs = ["macros.h"],
  deps = [
  ":bridge",
  ]
)

cc_library(
  name = "parse_command_line",
  hdrs = ["parse_command_line.h"],
  srcs = ["parse_command_line.cc"],
)

cc_library(
  name = "sequential_ht",
  hdrs = ["sequential_ht.h"],
  deps = [
  ":macros",
  "//pbbslib:utilities"
  ]
)

cc_library(
  name = "speculative_for",
  hdrs = ["speculative_for.h"],
  deps = [
  ":bridge"
  ]
)

cc_library(
  name = "union_find",
  hdrs = ["union_find.h"],
  srcs = ["union_find.cc"],
  deps = [
  ":macros",
  ":speculative_for"
  ]
)

cc_library(
  name = "vertex",
  hdrs = ["vertex.h"],
  deps = [
  ":macros",
  ":speculative_for"
  ]
)

cc_library(
  name = "vertex_subset",
  hdrs = ["vertex_subset.h"],
  deps = [
  ":bridge",
  ":macros",
  ]
)

cc_library(
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
