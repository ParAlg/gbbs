licenses(["notice"])
load("//internal_tools:build_defs.bzl", "gbbs_cc_test")
gbbs_cc_test(
    name = "regression_test",
    srcs = ["regression_test.cc"],
    deps = ["//benchmarks/PositiveWeightSSSP/DeltaStepping:DeltaStepping", "//gbbs:graph", "@googletest//:gtest_main"],
)
