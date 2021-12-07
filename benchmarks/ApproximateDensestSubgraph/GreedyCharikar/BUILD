licenses(["notice"])

cc_library(
    name = "DensestSubgraph",
    hdrs = ["DensestSubgraph.h"],
    deps = [
        "//benchmarks/KCore/JulienneDBS17:KCore",
        "//gbbs",
    ],
)

cc_binary(
    name = "DensestSubgraph_main",
    srcs = ["DensestSubgraph.cc"],
    deps = [":DensestSubgraph"],
)

package(
    default_visibility = ["//visibility:public"],
)
