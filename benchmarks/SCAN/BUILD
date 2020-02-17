cc_library(
    name = "scan",
    srcs = [
        "scan.cc",
        "scan_helpers.cc",
        "undirected_edge.cc",
    ],
    hdrs = [
        "scan.h",
        "scan_helpers.h",
        "undirected_edge.h",
    ],
    deps = [
        "//benchmarks/Connectivity/WorkEfficientSDB14:Connectivity",
        "//ligra:graph",
        "//ligra:macros",
        "//ligra/pbbslib:sparse_table",
        "//utils:assert",
    ],
)

cc_binary(
    name = "SCAN_main",
    srcs = ["run_scan.cc"],
    deps = [
        ":scan",
        "//ligra",
        "//utils:assert",
    ],
)

cc_test(
    name = "scan_test",
    srcs = ["unit_tests/scan_test.cc"],
    deps = [
        ":scan",
        "@googletest//:gtest_main",
    ],
)
