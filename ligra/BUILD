cc_library(
  name = "bridge",
  hdrs = ["bridge.h"],
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
  deps = [
  ":bridge",
  ":pbbs_strings",
  ]
)

cc_library(
  name = "graph_io",
  hdrs = ["graph_io.h"],
  deps = [
  ":bridge",
  ":graph",
  ":io",
  ":pbbs_strings",
  ]
)

cc_library(
  name = "bucket",
  hdrs = ["bucket.h"],
  deps = [
  ":maybe",
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
  deps = [
  ":bridge",
  ":edge_map_utils",
  ":graph",
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
  ":maybe",
  "//pbbslib:binary_search",
  "//pbbslib:utilities",
  ]
)

cc_library(
  name = "flags",
  hdrs = ["flags.h"],
)

cc_library(
  name = "graph",
  hdrs = ["graph.h"],
  deps = [
  ":bridge",
  ":compressed_vertex",
  ":flags",
  ":vertex",
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
  name = "histogram",
  hdrs = ["histogram.h"],
  deps = [
  ":bridge",
  ":sequential_ht",
  "//ligra/pbbslib:counting_sort_no_transpose",
  ]
)

cc_library(
  name = "ligra",
  hdrs = ["ligra.h"],
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

cc_library(
  name = "macros",
  hdrs = ["macros.h"],
  deps = [
  ":bridge",
  ]
)

cc_library(
  name = "maybe",
  hdrs = ["maybe.h"],
  deps = [
  ":macros",
  ]
)

cc_library(
  name = "parse_command_line",
  hdrs = ["parse_command_line.h"],
)

cc_library(
  name = "pbbs_strings",
  hdrs = ["pbbs_strings.h"],
  deps = ["//pbbslib/strings:string_basics"]
)

cc_library(
  name = "sequential_ht",
  hdrs = ["sequential_ht.h"],
  deps = [
  ":maybe",
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
  ":maybe"
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
