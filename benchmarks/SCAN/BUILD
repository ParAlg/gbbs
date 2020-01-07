cc_library(
    name = "scan",
    srcs = [
        "scan.cc",
        "undirected_edge.cc",
    ],
    hdrs = [
        "scan.h",
        "undirected_edge.h",
    ],
    deps = [
        "//ligra:graph",
        "//ligra:macros",
        "//ligra/pbbslib:sparse_table",
        "//utils:assert",
    ],
)

cc_binary(
    name = "run_scan",
    srcs = ["run_scan.cc"],
    deps = [
        ":scan",
        "//ligra",
        "//utils:assert",
    ],
)
