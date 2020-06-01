cc_library(
  name = "atomic_counter",
  hdrs = [
    "atomic_max_counter.h",
    "atomic_sum_counter.h",
  ],
  deps = [
    "//gbbs:macros",
    "//pbbslib:utilities",
  ],
)

cc_library(
  name = "counting_sort_no_transpose",
  hdrs = ["counting_sort_no_transpose.h"],
  deps = ["//gbbs:bridge"]
)

cc_library(
  name = "dyn_arr",
  hdrs = ["dyn_arr.h"],
  deps = ["//gbbs:bridge"]
)

cc_library(
  name = "resizable_table",
  hdrs = ["resizable_table.h"],
  deps = ["//gbbs:bridge"]
)

cc_library(
  name = "sparse_additive_map",
  hdrs = ["sparse_additive_map.h"],
  deps = ["//gbbs:bridge"]
)

cc_library(
  name = "sparse_table",
  hdrs = ["sparse_table.h"],
  deps = ["//gbbs:bridge"]
)

package(
  default_visibility = ["//visibility:public"],
)
