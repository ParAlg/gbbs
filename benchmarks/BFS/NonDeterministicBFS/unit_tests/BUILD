load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

gbbs_cc_test(
    name = "test_bfs",
    srcs = ["test_bfs.cc"],
    deps = [
        "//benchmarks/BFS/NonDeterministicBFS:BFS",
        "//ligra:graph",
        "//ligra:graph_test_utils",
        "//ligra:macros",
        "@googletest//:gtest_main",
    ],
)
