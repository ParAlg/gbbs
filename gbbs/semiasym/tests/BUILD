licenses(["notice"])

load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

gbbs_cc_test(
    name = "graph_filter_test",
    srcs = ["graph_filter_test.cc"],
    deps = [
        "//gbbs",
        "//gbbs:graph",
        "//gbbs/semiasym:graph_filter",
        "//gbbs/unit_tests:graph_test_utils",
        "@googletest//:gtest_main",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
