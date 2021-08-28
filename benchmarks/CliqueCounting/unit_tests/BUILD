load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

gbbs_cc_test(
    name = "test_clique",
    srcs = ["test_clique.cc"],
    deps = [
        "//benchmarks/CliqueCounting:Clique",
        "//gbbs:graph",
        "//gbbs:macros",
        "//gbbs/unit_tests:graph_test_utils",
        "@googletest//:gtest_main",
    ],
    size="small"
)
