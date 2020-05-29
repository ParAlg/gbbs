load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "scheduler",
  hdrs = ["scheduler.h"],
  srcs = ["scheduler.cc"],
  linkopts = ["-pthread"]
)

gbbs_cc_library(
  name = "parallel",
  hdrs = ["parallel.h"],
  srcs = ["parallel.cc"],
  deps = [":scheduler"]
)

gbbs_cc_library(
  name = "concurrent_stack",
  hdrs = ["concurrent_stack.h"],
)

gbbs_cc_library(
  name = "memory_size",
  hdrs = ["memory_size.h"],
  srcs = ["memory_size.cc"],
)

gbbs_cc_library(
  name = "alloc",
  hdrs = ["alloc.h"],
  srcs = ["alloc.cc"],
  deps = [
  ":concurrent_stack",
  ":memory_size",
  ":parallel",
  ],
)

gbbs_cc_library(
  name = "utilities",
  hdrs = ["utilities.h"],
  srcs = ["utilities.cc"],
  deps = [
  ":parallel",
  ":alloc"
  ],
)

gbbs_cc_library(
  name = "get_time",
  hdrs = ["get_time.h"],
  srcs = ["get_time.cc"],
)

gbbs_cc_library(
  name = "seq",
  hdrs = ["seq.h"],
  deps = [
  ":utilities"
  ]
)

gbbs_cc_library(
  name = "monoid",
  hdrs = ["monoid.h"],
)

gbbs_cc_library(
  name = "sequence_ops",
  hdrs = ["sequence_ops.h"],
  deps = [
  ":seq",
  ":monoid",
  ":utilities"
  ]
)

gbbs_cc_library(
  name = "random",
  hdrs = ["random.h"],
  deps = [
  ":utilities"
  ]
)

gbbs_cc_library(
  name = "transpose",
  hdrs = ["transpose.h"],
  deps = [
  ":utilities",
  ":get_time"
  ]
)

gbbs_cc_library(
  name = "counting_sort",
  hdrs = ["counting_sort.h"],
  deps = [
  ":utilities",
  ":sequence_ops",
  ":transpose"
  ]
)

gbbs_cc_library(
  name = "histogram",
  hdrs = ["histogram.h"],
  deps = [
  ":utilities",
  ":counting_sort"
  ]
)

gbbs_cc_library(
  name = "collect_reduce",
  hdrs = ["collect_reduce.h"],
  deps = [
  ":utilities",
  ":sequence_ops",
  ":transpose",
  ":histogram"
  ]
)

gbbs_cc_library(
  name = "quicksort",
  hdrs = ["quicksort.h"],
  deps = [
  ":utilities",
  ":sequence_ops"
  ]
)

gbbs_cc_library(
  name = "integer_sort",
  hdrs = ["integer_sort.h"],
  deps = [
  ":utilities",
  ":counting_sort",
  ":quicksort"
  ]
)

gbbs_cc_library(
  name = "binary_search",
  hdrs = ["binary_search.h"],
)

gbbs_cc_library(
  name = "merge",
  hdrs = ["merge.h"],
  deps = [
  ":utilities",
  ":seq",
  ":binary_search"
  ]
)

gbbs_cc_library(
  name = "merge_sort",
  hdrs = ["merge_sort.h"],
  deps = [
  ":utilities",
  ":quicksort",
  ":merge"
  ]
)

gbbs_cc_library(
  name = "bucket_sort",
  hdrs = ["bucket_sort.h"],
  deps = [
  ":utilities",
  ":sequence_ops",
  ":quicksort",
  ":get_time"
  ]
)

gbbs_cc_library(
  name = "sample_sort",
  hdrs = ["sample_sort.h"],
  deps = [
  ":utilities",
  ":sequence_ops",
  ":quicksort",
  ":merge_sort",
  ":transpose",
  ":bucket_sort",
  ":get_time"
  ]
)

gbbs_cc_library(
  name = "random_shuffle",
  hdrs = ["random_shuffle.h"],
  deps = [
  ":utilities",
  ":random",
  ":counting_sort"
  ]
)

gbbs_cc_library(
  name = "list_allocator",
  hdrs = ["list_allocator.h"],
  deps = [
  ":alloc",
  ":utilities",
  ":concurrent_stack",
  ":random_shuffle",
  ":memory_size"
  ]
)

gbbs_cc_library(
  name = "bag",
  hdrs = ["bag.h"],
  deps = [
  ":list_allocator"
  ]
)

gbbs_cc_library(
  name = "hash_table",
  hdrs = ["hash_table.h"],
  deps = [
  ":utilities",
  ":sequence_ops"
  ]
)

gbbs_cc_library(
  name = "sparse_mat_vec_mult",
  hdrs = ["sparse_mat_vec_mult.h"],
  deps = [
  ":utilities"
  ]
)

gbbs_cc_library(
  name = "kth_smallest",
  hdrs = ["kth_smallest.h"],
  deps = [
  ":sequence_ops",
  ":random"
  ]
)

gbbs_cc_library(
  name = "stlalgs",
  hdrs = ["stlalgs.h"],
  deps = [
  ":sequence_ops",
  ":sample_sort",
  ":kth_smallest"
  ]
)

gbbs_cc_library(
  name = "sequence",
  hdrs = ["sequence.h"],
  deps = [
  ":seq",
  ":sequence_ops",
  ":stlalgs"
  ]
)

gbbs_cc_library(
  name = "range_min",
  hdrs = ["range_min.h"],
  deps = [
  ":parallel",
  ":utilities",
  ":sequence"
  ]
)

gbbs_cc_library(
  name = "time_operations",
  hdrs = ["time_operations.h"],
  srcs = ["time_operations.cc"],
  deps = [
  ":get_time",
  ":random",
  ":counting_sort",
  ":collect_reduce",
  ":histogram",
  ":integer_sort",
  ":sample_sort",
  ":merge",
  ":merge_sort",
  ":bag",
  ":hash_table",
  ":sparse_mat_vec_mult",
  ":stlalgs",
  ":sequence_ops",
  ":monoid",
  ":range_min",
  ":utilities"
  ]
)

gbbs_cc_library(
  name = "assert",
  hdrs = ["assert.h"],
)

gbbs_cc_binary(
  name = "time_tests",
  srcs = ["time_tests.cc"],
  deps = [
   "//ligra:parse_command_line",
   ":alloc",
   ":get_time",
   ":sequence_ops",
   ":time_operations",
   ":utilities",
  ],
)

package(
  default_visibility = ["//visibility:public"],
)
