licenses(["notice"])

load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

gbbs_cc_test(
    name = "test_mis",
    srcs = ["test_mis.cc"],
    deps = [
        "//benchmarks/MaximalIndependentSet/RandomGreedy:MaximalIndependentSet",
        "//gbbs:graph",
        "//gbbs:macros",
        "//gbbs/unit_tests:graph_test_utils",
        "@googletest//:gtest_main",
    ],
)
