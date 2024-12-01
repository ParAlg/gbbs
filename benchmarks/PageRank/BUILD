licenses(["notice"])

load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "PageRank",
    hdrs = ["PageRank.h"],
    deps = [
        "//gbbs",
        "//gbbs:edge_map_reduce",
        "//gbbs/helpers:sparse_table",
    ],
)

cc_binary(
    name = "PageRank_main",
    srcs = ["PageRank.cc"],
    deps = [":PageRank"],
)

gbbs_cc_test(
    name = "PageRank_test",
    srcs = ["PageRank_test.cc"],
    deps = [
        "//benchmarks/PageRank:PageRank",
        "//gbbs:graph",
        "//gbbs:macros",
        "//gbbs/unit_tests:graph_test_utils",
        "@googletest//:gtest_main",
    ],
)
