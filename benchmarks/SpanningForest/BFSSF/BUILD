licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "SpanningForest",
    hdrs = ["SpanningForest.h"],
    deps = [
        "//benchmarks/SpanningForest:common",
        "//gbbs",
    ],
)

cc_binary(
    name = "SpanningForest_main",
    srcs = ["SpanningForest.cc"],
    deps = [
        ":SpanningForest",
        "//benchmarks/SpanningForest:common",
        "//gbbs",
    ],
)
