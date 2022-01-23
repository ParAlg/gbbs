licenses(["notice"])

cc_library(
    name = "MaximalIndependentSet",
    hdrs = ["MaximalIndependentSet.h"],
    deps = [
        "//gbbs",
        "//gbbs/helpers:speculative_for",
    ],
)

cc_binary(
    name = "MaximalIndependentSet_main",
    srcs = ["MaximalIndependentSet.cc"],
    deps = [":MaximalIndependentSet"],
)

package(
    default_visibility = ["//visibility:public"],
)
