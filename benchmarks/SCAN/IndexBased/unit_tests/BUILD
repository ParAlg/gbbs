load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

gbbs_cc_test(
    name = "scan_test",
    srcs = ["scan_test.cc"],
    deps = [
        "//benchmarks/SCAN/IndexBased:scan",
        "//ligra:graph",
        "//ligra:graph_test_utils",
        "//ligra:undirected_edge",
        "//ligra:vertex",
        "//pbbslib:seq",
        "@googletest//:gtest_main",
    ],
)

gbbs_cc_test(
    name = "utils_test",
    srcs = ["utils_test.cc"],
    deps = [
        "//benchmarks/SCAN/IndexBased:utils",
        "//ligra:graph_test_utils",
        "//pbbslib:seq",
        "@googletest//:gtest_main",
    ],
)
