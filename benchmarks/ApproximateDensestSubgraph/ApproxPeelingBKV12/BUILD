licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

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
