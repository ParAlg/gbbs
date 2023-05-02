licenses(["notice"])

load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

package(
    default_visibility = ["//visibility:public"],
)

gbbs_cc_test(
    name = "scan_test",
    srcs = ["scan_test.cc"],
    deps = [
        "//benchmarks/SCAN/Naive:scan",
        "//gbbs/helpers:undirected_edge",
        "//gbbs/unit_tests:graph_test_utils",
        "@googletest//:gtest_main",
    ],
)
