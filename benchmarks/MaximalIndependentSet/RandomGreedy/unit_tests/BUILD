load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

gbbs_cc_test(
    name = "test_mis",
    srcs = ["test_mis.cc"],
    deps = [
        "//benchmarks/MaximalIndependentSet/RandomGreedy:MaximalIndependentSet",
        "//gbbs:graph",
        "//gbbs/unit_tests:graph_test_utils",
        "//gbbs:macros",
        "@googletest//:gtest_main",
    ],
)
