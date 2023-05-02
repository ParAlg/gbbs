licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "Spanner",
    hdrs = ["Spanner.h"],
    deps = [
        "//benchmarks/LowDiameterDecomposition/MPX13:LowDiameterDecomposition",
        "//gbbs",
        "//gbbs/helpers:sparse_table",
    ],
)

cc_binary(
    name = "Spanner_main",
    srcs = ["Spanner.cc"],
    deps = [":Spanner"],
)
