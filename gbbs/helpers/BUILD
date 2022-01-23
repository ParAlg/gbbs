licenses(["notice"])

cc_library(
    name = "assert",
    hdrs = ["assert.h"],
    deps = [],
)

cc_library(
    name = "atomic_counter",
    hdrs = [
        "atomic_max_counter.h",
        "atomic_sum_counter.h",
    ],
    deps = [
        "//gbbs:bridge",
        "//gbbs:macros",
    ],
)

cc_library(
    name = "counting_sort_no_transpose",
    hdrs = ["counting_sort_no_transpose.h"],
    deps = ["//gbbs:bridge"],
)

cc_library(
    name = "dyn_arr",
    hdrs = ["dyn_arr.h"],
    deps = ["//gbbs:bridge"],
)

cc_library(
    name = "histogram",
    hdrs = ["histogram.h"],
    deps = [
        ":sequential_ht",
        "//gbbs:bridge",
        "//gbbs:macros",
        "//gbbs/helpers:counting_sort_no_transpose",
    ],
)

cc_library(
    name = "parse_command_line",
    hdrs = ["parse_command_line.h"],
)

cc_library(
    name = "resizable_table",
    hdrs = ["resizable_table.h"],
    deps = ["//gbbs:bridge"],
)

cc_library(
    name = "sequential_ht",
    hdrs = ["sequential_ht.h"],
    deps = [
        "//gbbs:bridge",
        "//gbbs:macros",
    ],
)

cc_library(
    name = "sparse_additive_map",
    hdrs = ["sparse_additive_map.h"],
    deps = ["//gbbs:bridge"],
)

cc_library(
    name = "sparse_table",
    hdrs = ["sparse_table.h"],
    deps = ["//gbbs:bridge"],
)

cc_library(
    name = "speculative_for",
    hdrs = ["speculative_for.h"],
    deps = [
        "//gbbs:bridge",
    ],
)

cc_library(
    name = "undirected_edge",
    srcs = ["undirected_edge.cc"],
    hdrs = ["undirected_edge.h"],
    deps = [
        "//gbbs:macros",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
