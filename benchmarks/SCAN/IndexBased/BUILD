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
        "//benchmarks/Connectivity/UnionFind:union_find_rules",
        "//benchmarks/TriangleCounting/ShunTangwongsan15:Triangle_main",
        "//ligra:graph",
        "//ligra:macros",
        "//ligra:undirected_edge",
        "//ligra/pbbslib:sparse_table",
        "//pbbslib:get_time",
        "//pbbslib:parallel",
        "//pbbslib:utilities",
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
