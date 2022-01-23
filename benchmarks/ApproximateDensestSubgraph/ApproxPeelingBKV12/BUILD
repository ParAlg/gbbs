licenses(["notice"])

cc_library(
    name = "DensestSubgraph",
    hdrs = ["DensestSubgraph.h"],
    deps = [
        "//gbbs",
        "//gbbs:edge_map_reduce",
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
