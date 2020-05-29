cc_library(
    name = "scan",
    srcs = [
        "scan.cc",
        "scan_helpers.cc",
    ],
    hdrs = [
        "intersect.h",
        "scan.h",
        "scan_helpers.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":utils",
        "//benchmarks/Connectivity/UnionFind:union_find_rules",
        "//benchmarks/TriangleCounting/ShunTangwongsan15:Triangle_main",
        "//ligra:bridge",
        "//ligra:graph",
        "//ligra:graph_mutation",
        "//ligra:macros",
        "//ligra:undirected_edge",
        "//ligra/pbbslib:sparse_table",
        "//pbbslib:binary_search",
        "//pbbslib:get_time",
        "//pbbslib:parallel",
        "//pbbslib:sample_sort",
        "//pbbslib:seq",
        "//pbbslib:utilities",
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
        "//ligra:bridge",
        "//ligra:graph",
        "//ligra:macros",
        "//pbbslib:collect_reduce",
        "//pbbslib:monoid",
        "//pbbslib:seq",
    ],
)

cc_binary(
    name = "SCAN_main",
    srcs = ["run_scan.cc"],
    deps = [
        ":scan",
        "//ligra",
    ],
)
