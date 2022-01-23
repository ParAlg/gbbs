licenses(["notice"])

package(default_visibility = ["//visibility:public"])

load("//internal_tools:py_extension.bzl", "py_extension")

cc_library(
    name = "BFS_lib",
    srcs = ["BFS_lib.cc"],
    hdrs = ["BFS_lib.h"],
    deps = [
        "//benchmarks/BFS/NonDeterministicBFS:BFS",
        "//gbbs",
    ],
)

cc_library(
    name = "CC_lib",
    srcs = ["CC_lib.cc"],
    hdrs = ["CC_lib.h"],
    deps = [
        "//benchmarks/Connectivity/WorkEfficientSDB14:Connectivity",
        "//gbbs",
    ],
)

cc_library(
    name = "KCore_lib",
    srcs = ["KCore_lib.cc"],
    hdrs = ["KCore_lib.h"],
    deps = [
        "//benchmarks/KCore/JulienneDBS17:KCore",
        "//gbbs",
    ],
)

cc_library(
    name = "HAC_lib",
    srcs = ["HAC_lib.cc"],
    hdrs = ["HAC_lib.h"],
    deps = [
        "//benchmarks/Clustering/SeqHAC:HAC_api",
        "//gbbs",
    ],
)

cc_library(
    name = "BellmanFord_lib",
    srcs = ["BellmanFord_lib.cc"],
    hdrs = ["BellmanFord_lib.h"],
    deps = [
        "//benchmarks/GeneralWeightSSSP/BellmanFord",
        "//gbbs",
    ],
)

cc_library(
    name = "DeltaStepping_lib",
    srcs = ["DeltaStepping_lib.cc"],
    hdrs = ["DeltaStepping_lib.h"],
    deps = [
        "//benchmarks/PositiveWeightSSSP/DeltaStepping",
        "//gbbs",
    ],
)

cc_library(
    name = "PageRank_lib",
    srcs = ["PageRank_lib.cc"],
    hdrs = ["PageRank_lib.h"],
    deps = [
        "//benchmarks/PageRank",
        "//gbbs",
    ],
)

cc_library(
    name = "MinimumSpanningForest_lib",
    srcs = ["MinimumSpanningForest_lib.cc"],
    hdrs = ["MinimumSpanningForest_lib.h"],
    deps = [
        "//benchmarks/MinimumSpanningForest/Boruvka:MinimumSpanningForest",
        #    "//benchmarks/MinimumSpanningForest/Kruskal:MinimumSpanningForest",
        "//gbbs",
    ],
)

cc_library(
    name = "ApproximateSetCover_lib",
    srcs = ["ApproximateSetCover_lib.cc"],
    hdrs = ["ApproximateSetCover_lib.h"],
    deps = [
        "//benchmarks/ApproximateSetCover/MANISBPT11:ApproximateSetCover",
        "//gbbs",
    ],
)

py_extension(
    name = "gbbs_lib",
    srcs = ["gbbs_lib.cc"],
    deps = [
        ":ApproximateSetCover_lib",
        ":BFS_lib",
        ":BellmanFord_lib",
        ":CC_lib",
        ":DeltaStepping_lib",
        ":HAC_lib",
        ":KCore_lib",
        ":MinimumSpanningForest_lib",
        ":PageRank_lib",
        "//gbbs",
        "//gbbs:compressed_vertex",
        "//gbbs:graph",
        "//gbbs:graph_io",
        "//gbbs:vertex",
        "//gbbs:vertex_subset",
        "@local_config_python//:python_headers",
        "@pybind11",
    ],
)
