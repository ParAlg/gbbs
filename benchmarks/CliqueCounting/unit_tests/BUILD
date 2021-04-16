load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

gbbs_cc_test(
    name = "test_clique",
    srcs = ["test_clique.cc"],
    deps = [
        "//benchmarks/CliqueCounting:Clique",
        "//gbbs:graph",
        "//gbbs:graph_test_utils",
        "//gbbs:macros",
        "@googletest//:gtest_main",
    ],
    size="small"
)
