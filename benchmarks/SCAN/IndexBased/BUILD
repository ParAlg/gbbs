licenses(["notice"])

cc_library(
    name = "scan",
    srcs = [
        "scan.cc",
        "scan_helpers.cc",
    ],
    hdrs = [
        "scan.h",
        "scan_helpers.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":similarity_measure",
        ":utils",
        "//benchmarks/Connectivity/UnionFind:union_find_rules",
        "//gbbs:bridge",
        "//gbbs:graph",
        "//gbbs:macros",
        "//gbbs/helpers:sparse_table",
    ],
)

cc_library(
    name = "similarity_measure",
    srcs = [
        "similarity_measure.cc",
    ],
    hdrs = [
        "intersect.h",
        "similarity_measure.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//benchmarks/TriangleCounting/ShunTangwongsan15:Triangle",
        "//gbbs:bridge",
        "//gbbs:graph_mutation",
        "//gbbs:macros",
        "//gbbs:vertex",
        "//gbbs/helpers:assert",
    ],
)

cc_library(
    name = "utils",
    srcs = [
        "utils.cc",
    ],
    hdrs = [
        "utils.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//gbbs:bridge",
        "//gbbs:graph",
        "//gbbs:macros",
    ],
)

cc_binary(
    name = "SCAN_main",
    srcs = ["run_scan.cc"],
    deps = [
        ":scan",
        ":similarity_measure",
        ":utils",
        "//gbbs",
    ],
)
