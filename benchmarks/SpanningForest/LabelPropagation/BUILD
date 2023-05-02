licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "SpanningForest",
    hdrs = ["SpanningForest.h"],
    deps = [
        "//benchmarks/Connectivity:common",
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
        "//benchmarks/SpanningForest/BFSSF:SpanningForest",
        "//gbbs",
    ],
)
