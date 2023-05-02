licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "PageRank",
    hdrs = ["PageRank.h"],
    deps = [
        "//gbbs",
        "//gbbs:edge_map_reduce",
        "//gbbs/helpers:sparse_table",
    ],
)

cc_binary(
    name = "PageRank_main",
    srcs = ["PageRank.cc"],
    deps = [":PageRank"],
)
