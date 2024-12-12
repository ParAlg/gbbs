licenses(["notice"])

load("//internal_tools:build_defs.bzl", "gbbs_cc_test")

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "LabelPropagation",
    hdrs = ["LabelPropagation.h"],
    deps = [
        "//gbbs",
        "//benchmarks/GraphColoring/Hasenplaugh14:GraphColoring"
    ],
)

cc_binary(
    name = "LabelPropagation_main",
    srcs = ["LabelPropagation.cc"],
    deps = [":LabelPropagation"],
)

#gbbs_cc_test(
#    name = "LabelPropagation_test",
#    srcs = ["LabelPropagation_test.cc"],
#    deps = [
#        "//benchmarks/LabelPropagation:LabelPropagation",
#        "//gbbs:graph",
#        "//gbbs:macros",
#        "//gbbs/unit_tests:graph_test_utils",
#        "@googletest//:gtest_main",
#    ],
#)
