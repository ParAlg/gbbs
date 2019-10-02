cc_library(
  name = "scheduler",
  hdrs = ["scheduler.h"],
  linkopts = ["-pthread"]
)

cc_library(
  name = "parallel",
  hdrs = ["parallel.h"],
  deps = [":scheduler"]
)

cc_library(
  name = "concurrent_stack",
  hdrs = ["concurrent_stack.h"],
)

cc_library(
  name = "memory_size",
  hdrs = ["memory_size.h"],
)

cc_library(
  name = "alloc",
  hdrs = ["alloc.h"],
  deps = [
  ":concurrent_stack",
  ":memory_size",
  ],
)

cc_library(
  name = "utilities",
  hdrs = ["utilities.h"],
  deps = [
  ":parallel",
  ":alloc"
  ],
)

cc_library(
  name = "block_allocator",
  hdrs = ["block_allocator.h"],
  deps = [
  ":concurrent_stack",
  ":memory_size",
  ":utilities"
  ],
)

cc_library(
  name = "allocator",
  hdrs = ["allocator.h"],
  deps = [
  ":block_allocator",
  ],
)

cc_library(
  name = "get_time",
  hdrs = ["get_time.h"],
)

cc_library(
  name = "parse_command_line",
  hdrs = ["parse_command_line.h"],
)


cc_library(
  name = "seq",
  hdrs = ["seq.h"],
  deps = [
  ":utilities"
  ]
)

cc_library(
  name = "monoid",
  hdrs = ["monoid.h"],
)

cc_library(
  name = "sequence_ops",
  hdrs = ["sequence_ops.h"],
  deps = [
  ":seq",
  ":monoid",
  ":utilities"
  ]
)

cc_test(
  name = "sequence_test",
  srcs = ["sequence_test.cpp"],
  copts = ["-Iexternal/gtest/include"],
  deps = [
      "@googletest//:gtest_main",
  ],
)

cc_library(
  name = "random",
  hdrs = ["random.h"],
  deps = [
  ":utilities"
  ]
)

cc_library(
  name = "transpose",
  hdrs = ["transpose.h"],
  deps = [
  ":utilities",
  ":get_time"
  ]
)

cc_library(
  name = "counting_sort",
  hdrs = ["counting_sort.h"],
  deps = [
  ":utilities",
  ":sequence_ops",
  ":transpose"
  ]
)

cc_library(
  name = "histogram",
  hdrs = ["histogram.h"],
  deps = [
  ":utilities",
  ":counting_sort"
  ]
)

cc_library(
  name = "collect_reduce",
  hdrs = ["collect_reduce.h"],
  deps = [
  ":utilities",
  ":sequence_ops",
  ":transpose",
  ":histogram"
  ]
)

cc_library(
  name = "quicksort",
  hdrs = ["quicksort.h"],
  deps = [
  ":utilities",
  ":sequence_ops"
  ]
)

cc_library(
  name = "integer_sort",
  hdrs = ["integer_sort.h"],
  deps = [
  ":utilities",
  ":counting_sort",
  ":quicksort"
  ]
)

cc_library(
  name = "binary_search",
  hdrs = ["binary_search.h"],
)

cc_library(
  name = "merge",
  hdrs = ["merge.h"],
  deps = [
  ":utilities",
  ":seq",
  ":binary_search"
  ]
)

cc_library(
  name = "merge_sort",
  hdrs = ["merge_sort.h"],
  deps = [
  ":utilities",
  ":quicksort",
  ":merge"
  ]
)

cc_library(
  name = "bucket_sort",
  hdrs = ["bucket_sort.h"],
  deps = [
  ":utilities",
  ":sequence_ops",
  ":quicksort",
  ":get_time"
  ]
)

cc_library(
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

cc_library(
  name = "random_shuffle",
  hdrs = ["random_shuffle.h"],
  deps = [
  ":utilities",
  ":random",
  ":counting_sort"
  ]
)

cc_library(
  name = "list_allocator",
  hdrs = ["list_allocator.h"],
  deps = [
  ":utilities",
  ":concurrent_stack",
  ":random_shuffle",
  ":memory_size"
  ]
)

cc_library(
  name = "bag",
  hdrs = ["bag.h"],
  deps = [
  ":list_allocator"
  ]
)

cc_library(
  name = "hash_table",
  hdrs = ["hash_table.h"],
  deps = [
  ":utilities",
  ":sequence_ops"
  ]
)

cc_library(
  name = "sparse_mat_vec_mult",
  hdrs = ["sparse_mat_vec_mult.h"],
  deps = [
  ":utilities"
  ]
)

cc_library(
  name = "kth_smallest",
  hdrs = ["kth_smallest.h"],
  deps = [
  ":sequence_ops",
  ":random"
  ]
)

cc_library(
  name = "stlalgs",
  hdrs = ["stlalgs.h"],
  deps = [
  ":sequence_ops",
  ":sample_sort",
  ":kth_smallest"
  ]
)

cc_library(
  name = "sequence",
  hdrs = ["sequence.h"],
  deps = [
  ":seq",
  ":sequence_ops",
  ":stlalgs"
  ]
)

cc_library(
  name = "range_min",
  hdrs = ["range_min.h"],
  deps = [
  ":parallel",
  ":utilities",
  ":sequence"
  ]
)

cc_library(
  name = "time_operations",
  hdrs = ["time_operations.h"],
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


cc_binary(
  name = "time_tests",
  srcs = ["time_tests.cc"],
  deps = [
   ":alloc",
   ":allocator",
   ":get_time",
   ":parse_command_line",
   ":sequence_ops",
   ":time_operations",
   ":utilities",
  ],
)

package(
  default_visibility = ["//visibility:public"],
)
