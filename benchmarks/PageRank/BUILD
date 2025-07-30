licenses(["notice"])

load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "PageRank",
    hdrs = ["PageRank.h"],
    deps = [
        "//gbbs:bridge",
        "//gbbs:edge_map_data",
        "//gbbs:edge_map_reduce",
        "//gbbs:flags",
        "//gbbs:macros",
        "//gbbs:vertex_subset",
        "//gbbs/helpers:assert",
        "@parlaylib//parlay:monoid",
        "@parlaylib//parlay:sequence",
    ],
)

cc_binary(
    name = "PageRank_main",
    srcs = ["PageRank.cc"],
    deps = [
        ":PageRank",
        "//gbbs:benchmark",
        "//gbbs:bridge",
        "//gbbs/helpers:parse_command_line",
    ],
)

gbbs_cc_test(
    name = "PageRank_test",
    srcs = ["PageRank_test.cc"],
    deps = [
        ":PageRank",
        "//gbbs:bridge",
        "//gbbs:graph",
        "//gbbs:macros",
        "//gbbs/helpers:undirected_edge",
        "//gbbs/unit_tests:graph_test_utils",
        "@googletest//:gtest_main",
    ],
    size = "small",
)