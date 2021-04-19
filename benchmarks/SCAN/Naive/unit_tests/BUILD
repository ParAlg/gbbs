load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

gbbs_cc_test(
    name = "scan_test",
    srcs = ["scan_test.cc"],
    deps = [
        "//benchmarks/SCAN/Naive:scan",
        "//gbbs:graph_test_utils",
        "//gbbs:undirected_edge",
        "@googletest//:gtest_main",
    ],
)
