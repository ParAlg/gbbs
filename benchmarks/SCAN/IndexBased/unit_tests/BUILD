load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

gbbs_cc_test(
    name = "scan_test",
    srcs = ["scan_test.cc"],
    deps = [
        ":similarity_measure_test_utils",
        "//benchmarks/SCAN/IndexBased:scan",
        "//benchmarks/SCAN/IndexBased:similarity_measure",
        "//gbbs:bridge",
        "//gbbs:graph",
        "//gbbs:graph_test_utils",
        "//gbbs:undirected_edge",
        "//gbbs:vertex",
        "@googletest//:gtest_main",
    ],
)

#gbbs_cc_test(
#    name = "similarity_measure_test",
#    srcs = ["similarity_measure_test.cc"],
#    deps = [
#        ":similarity_measure_test_utils",
#        "//benchmarks/SCAN/IndexBased:similarity_measure",
#        "//gbbs:graph_test_utils",
#        "//gbbs:macros",
#        "@googletest//:gtest_main",
#    ],
#)
#
#gbbs_cc_test(
#    name = "utils_test",
#    srcs = ["utils_test.cc"],
#    deps = [
#        "//benchmarks/SCAN/IndexBased:utils",
#        "//gbbs:graph_test_utils",
#        "@googletest//:gtest_main",
#    ],
#)

cc_library(
    name = "similarity_measure_test_utils",
    hdrs = ["similarity_measure_test_utils.h"],
    deps = [
        "//benchmarks/SCAN/IndexBased:similarity_measure",
        "@googletest//:gtest",
    ],
)
