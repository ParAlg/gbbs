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
        "@abseil-cpp//absl/container:flat_hash_set",
        "@abseil-cpp//absl/types:span",
        "@parlaylib//parlay:monoid",
        "@parlaylib//parlay:sequence",
    ],
)

cc_binary(
    name = "PageRank_main",
    srcs = ["PageRank.cc"],
    deps = [
        ":PageRank",
        ":PageRank_delta",
        ":PageRank_edgeMapReduce",
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
        ":PageRank_edgeMapReduce",
        "//gbbs:bridge",
        "//gbbs:graph",
        "//gbbs:macros",
        "//gbbs/helpers:undirected_edge",
        "//gbbs/unit_tests:graph_test_utils",
        "@googletest//:gtest_main",
    ],
    size = "small",
)

cc_library(
    name = "PageRank_delta",
    hdrs = ["PageRank_delta.h"],
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

cc_library(
    name = "PageRank_edgeMapReduce",
    hdrs = ["PageRank_edgeMapReduce.h"],
    deps = [
        "//gbbs:bridge",
        "//gbbs:edge_map_reduce",
        "//gbbs:flags",
        "//gbbs:macros",
        "//gbbs:vertex_subset",
        "//gbbs/helpers:assert",
        "@abseil-cpp//absl/container:flat_hash_set",
        "@parlaylib//parlay:monoid",
        "@parlaylib//parlay:sequence",
    ],
)