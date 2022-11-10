licenses(["notice"])

load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

gbbs_cc_test(
    name = "test_max_flow",
    srcs = ["test_max_flow.cc"],
    deps = [
        "//benchmarks/MaximumFlow/PushRelabelBSB15:PushRelabel",
        "//gbbs:graph",
        "//gbbs:macros",
        "//gbbs/unit_tests:graph_test_utils",
        "@googletest//:gtest_main",
    ],
)
