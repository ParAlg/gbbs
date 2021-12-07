licenses(["notice"])

load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

gbbs_cc_test(
    name = "test_cc",
    srcs = ["test_cc.cc"],
    deps = [
        "//benchmarks/Connectivity/WorkEfficientSDB14:Connectivity",
        "//gbbs:graph",
        "//gbbs:macros",
        "//gbbs/unit_tests:graph_test_utils",
        "@googletest//:gtest_main",
    ],
)
