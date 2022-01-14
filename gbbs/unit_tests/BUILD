licenses(["notice"])

load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

package(
    default_visibility = ["//visibility:public"],
)

gbbs_cc_test(
    name = "graph_io_test",
    srcs = ["graph_io_test.cc"],
    deps = [
        ":graph_test_utils",
        "//gbbs:graph_io",
        "@googletest//:gtest_main",
    ],
)

gbbs_cc_test(
    name = "graph_test",
    srcs = ["graph_test.cc"],
    deps = [
        "//gbbs:graph",
        "@googletest//:gtest_main",
        "@parlaylib//parlay:primitives",
    ],
)

gbbs_cc_test(
    name = "undirected_edge_test",
    srcs = ["undirected_edge_test.cc"],
    deps = [
        "//gbbs/helpers:undirected_edge",
        "@googletest//:gtest_main",
    ],
)

cc_library(
    name = "graph_test_utils",
    testonly = 1,
    srcs = ["graph_test_utils.cc"],
    hdrs = ["graph_test_utils.h"],
    deps = [
        "//gbbs:bridge",
        "//gbbs:graph",
        "//gbbs:macros",
        "//gbbs:vertex",
        "//gbbs/helpers:undirected_edge",
        "@googletest//:gtest",
    ],
)
