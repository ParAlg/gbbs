load("//internal_tools:build_defs.bzl", "gbbs_cc_library")

gbbs_cc_library(
  name = "atomic_counter",
  hdrs = [
    "atomic_max_counter.h",
    "atomic_sum_counter.h",
  ],
  deps = [
    "//ligra:macros",
    "//pbbslib:utilities",
  ],
)

gbbs_cc_library(
  name = "counting_sort_no_transpose",
  hdrs = ["counting_sort_no_transpose.h"],
  deps = ["//ligra:bridge"]
)

gbbs_cc_library(
  name = "dyn_arr",
  hdrs = ["dyn_arr.h"],
  deps = ["//ligra:bridge"]
)

gbbs_cc_library(
  name = "resizable_table",
  hdrs = ["resizable_table.h"],
  deps = ["//ligra:bridge"]
)

gbbs_cc_library(
  name = "sparse_additive_map",
  hdrs = ["sparse_additive_map.h"],
  deps = ["//ligra:bridge"]
)

gbbs_cc_library(
  name = "sparse_table",
  hdrs = ["sparse_table.h"],
  deps = ["//ligra:bridge"]
)

package(
  default_visibility = ["//visibility:public"],
)
