licenses(["notice"])

load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

gbbs_cc_test(
    name = "test_maximal_matching",
    srcs = ["test_maximal_matching.cc"],
    deps = [
        "//benchmarks/MaximalMatching/RandomGreedy:MaximalMatching",
        "//gbbs:graph",
        "//gbbs/unit_tests:graph_test_utils",
        "//gbbs:macros",
        "@googletest//:gtest_main",
    ],
)
