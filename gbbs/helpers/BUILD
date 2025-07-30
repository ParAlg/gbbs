licenses(["notice"])

load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

package(
    default_visibility = ["//visibility:public"],
)

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
    name = "directed_edge",
    srcs = ["directed_edge.cc"],
    hdrs = ["directed_edge.h"],
    deps = [
        "//gbbs:macros",
    ],
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
    name = "parallel_for_with_status",
    hdrs = ["parallel_for_with_status.h"],
    deps = [
        "//gbbs:bridge",
        "@abseil-cpp//absl/status",
        "@parlaylib//parlay:parallel",
    ],
)

gbbs_cc_test(
    name = "parallel_for_with_status_test",
    srcs = ["parallel_for_with_status_test.cc"],
    deps = [
        ":parallel_for_with_status",
        ":status_macros",
        "@googletest//:gtest_main",
        "@abseil-cpp//absl/base:core_headers",
        "@abseil-cpp//absl/log",
        "@abseil-cpp//absl/log:absl_log",
        "@abseil-cpp//absl/log:check",
        "@abseil-cpp//absl/status",
        "@abseil-cpp//absl/strings",
        "@abseil-cpp//absl/synchronization",
        "@google_benchmark//:benchmark",
        "@parlaylib//parlay:parallel",
    ],
)

cc_library(
    name = "parse_command_line",
    hdrs = ["parse_command_line.h"],
)

cc_library(
    name = "progress_reporting",
    srcs = ["progress_reporting.cc"],
    hdrs = ["progress_reporting.h"],
    deps = [
        "@abseil-cpp//absl/functional:any_invocable",
        "@abseil-cpp//absl/log:absl_check",
        "@abseil-cpp//absl/status",
        "@abseil-cpp//absl/strings",
    ],
)

gbbs_cc_test(
    name = "progress_reporting_test",
    srcs = ["progress_reporting_test.cc"],
    deps = [
        ":progress_reporting",
        ":progress_reporting_mock",
        ":status_macros",
        "@abseil-cpp//absl/status",
        "@googletest//:gtest_main",
    ],
)

cc_library(
    name = "progress_reporting_mock",
    testonly = True,
    hdrs = ["progress_reporting_mock.h"],
    deps = [
        ":progress_reporting",
        "@abseil-cpp//absl/base:core_headers",
    ],
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
    name = "status_macros",
    srcs = ["status_macros.cc"],
    hdrs = ["status_macros.h"],
    deps = [
        "@googletest//:gtest",
        "@abseil-cpp//absl/log:absl_log",
        "@abseil-cpp//absl/status",
        "@abseil-cpp//absl/status:statusor",
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
