load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "GraphColoring",
    hdrs = ["GraphColoring.h"],
    deps = [
        "//gbbs:bridge",
        "//gbbs:edge_map_data",
        "//gbbs:flags",
        "//gbbs:macros",
        "//gbbs:vertex_subset",
        "@abseil-cpp//absl/types:span",
        "@parlaylib//parlay:random",
    ],
)

cc_binary(
    name = "GraphColoring_main",
    srcs = ["GraphColoring.cc"],
    deps = [
        ":GraphColoring",
        "//gbbs:benchmark",
        "//gbbs/helpers:parse_command_line",
    ],
)

gbbs_cc_test(
    name = "GraphColoring_test",
    srcs = ["GraphColoring_test.cc"],
    deps = [
        ":GraphColoring",
        "//gbbs:bridge",
        "//gbbs:graph",
        "//gbbs:graph_generator",
        "//gbbs:macros",
        "//gbbs:vertex",
        "@googletest//:gtest_main",
        "@parlaylib//parlay:monoid",
        "@parlaylib//parlay:primitives",
    ],
)

gbbs_cc_test(
    name = "GraphColoring_fuzz_test",
    srcs = ["GraphColoring_fuzz_test.cc"],
    tags = ["fuzz_test"],
    deps = [
        ":GraphColoring",
        "//gbbs:bridge",
        "//gbbs:graph",
        "//gbbs:macros",
        "//gbbs:vertex",
        "@fuzztest//fuzztest",
        "@fuzztest//fuzztest:fuzztest_gtest_main",
        "@googletest//:gtest",
    ],
)
