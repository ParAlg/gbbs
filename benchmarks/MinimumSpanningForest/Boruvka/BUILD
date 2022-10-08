licenses(["notice"])

cc_library(
    name = "MinimumSpanningForest",
    hdrs = ["MinimumSpanningForest.h"],
    deps = [
        "//benchmarks/SpanningForest/SDB14:SpanningForest",
        "//gbbs",
        "//gbbs/helpers:sparse_table",
    ],
)

cc_binary(
    name = "MinimumSpanningForest_main",
    srcs = ["MinimumSpanningForest.cc"],
    deps = [":MinimumSpanningForest"],
)

package(
    default_visibility = ["//visibility:public"],
)
