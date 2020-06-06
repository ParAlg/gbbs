load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

gbbs_cc_test(
    name = "seq_test",
    srcs = ["seq_test.cc"],
    deps = [
        "//pbbslib:seq",
        "@googletest//:gtest_main",
    ],
)

gbbs_cc_test(
    name = "sequence_ops_test",
    srcs = ["sequence_ops_test.cc"],
    deps = [
        "//pbbslib:monoid",
        "//pbbslib:seq",
        "//pbbslib:sequence_ops",
        "@googletest//:gtest_main",
    ],
)
