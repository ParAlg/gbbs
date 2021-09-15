licenses(["notice"])

cc_library(
    name = "HAC_api",
    hdrs = ["HAC_api.h"],
    deps = [
        ":HAC_configuration",
        ":HeapBased",
        ":NNChainBased",
        "//benchmarks/Clustering/SeqHAC/AvgLinkageUtils:HeapBased",
        "//gbbs",
    ],
)

cc_library(
    name = "ClusteredGraph",
    hdrs = ["ClusteredGraph.h"],
    deps = [
        "//gbbs",
        "@PAM//pam",
    ],
)

cc_library(
    name = "HAC_configuration",
    hdrs = ["HAC_configuration.h"],
    deps = [
        "//gbbs:macros",
    ],
)

cc_library(
    name = "HeapBased",
    hdrs = ["HeapBased.h"],
    deps = [
        ":ClusteredGraph",
        "//gbbs:macros",
    ],
)

cc_library(
    name = "NNChainBased",
    hdrs = ["NNChainBased.h"],
    deps = [
        ":ClusteredGraph",
        "//gbbs:macros",
    ],
)

cc_binary(
    name = "HACDissimilarity",
    srcs = ["HACDissimilarity.cc"],
    deps = [
        ":HAC_configuration",
        ":HeapBased",
        ":NNChainBased",
        "//benchmarks/Clustering/SeqHAC/AvgLinkageUtils:HeapBased",
    ],
)

cc_binary(
    name = "HACSimilarity",
    srcs = ["HACSimilarity.cc"],
    deps = [
        ":HAC_configuration",
        ":HeapBased",
        ":NNChainBased",
        "//benchmarks/Clustering/SeqHAC/AvgLinkageUtils:HeapBased",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
