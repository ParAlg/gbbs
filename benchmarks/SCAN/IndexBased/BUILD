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
        "//benchmarks/Connectivity/WorkEfficientSDB14:Connectivity",
        "//ligra:graph",
        "//ligra:graph_test_utils",
        "//ligra:macros",
        "//ligra:undirected_edge",
        "//ligra/pbbslib:sparse_table",
        "//pbbslib:get_time",
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
