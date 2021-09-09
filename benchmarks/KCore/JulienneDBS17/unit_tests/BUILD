licenses(["notice"])

load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

gbbs_cc_test(
    name = "test_kcore",
    srcs = ["test_kcore.cc"],
    deps = [
        "//benchmarks/KCore/JulienneDBS17:KCore",
        "//gbbs:graph",
        "//gbbs/unit_tests:graph_test_utils",
        "//gbbs:macros",
        "@googletest//:gtest_main",
    ],
)
