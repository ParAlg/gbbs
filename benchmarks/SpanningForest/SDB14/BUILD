licenses(["notice"])

cc_library(
    name = "contract_sf",
    hdrs = ["contract_sf.h"],
    deps = [
        "//gbbs:graph",
        "//gbbs:macros",
        "//gbbs/helpers:sparse_table",
    ],
)

cc_library(
    name = "SpanningForest",
    hdrs = ["SpanningForest.h"],
    deps = [
        ":contract_sf",
        "//benchmarks/LowDiameterDecomposition/MPX13:LowDiameterDecomposition",
        "//gbbs",
        "//gbbs/helpers:sparse_table",
    ],
)

cc_binary(
    name = "SpanningForest_main",
    srcs = ["SpanningForest.cc"],
    deps = [
        ":SpanningForest",
        "//benchmarks/SpanningForest:common",
        "//benchmarks/SpanningForest/BFSSF:SpanningForest",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
