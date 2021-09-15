licenses(["notice"])

load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

gbbs_cc_test(
    name = "test_bfs",
    srcs = ["test_bfs.cc"],
    deps = [
        "//benchmarks/BFS/NonDeterministicBFS:BFS",
        "//gbbs:graph",
        "//gbbs:macros",
        "//gbbs/unit_tests:graph_test_utils",
        "@googletest//:gtest_main",
    ],
)
