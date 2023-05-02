licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "MinimumSpanningForest",
    hdrs = ["MinimumSpanningForest.h"],
    deps = [
        "//benchmarks/Connectivity:common",
        "//benchmarks/Connectivity/UnionFind:union_find_rules",
        "//gbbs",
        "//gbbs/helpers:sparse_table",
    ],
)

cc_binary(
    name = "MinimumSpanningForest_main",
    srcs = ["MinimumSpanningForest.cc"],
    deps = [":MinimumSpanningForest"],
)
