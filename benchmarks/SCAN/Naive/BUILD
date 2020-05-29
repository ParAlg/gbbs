load("//internal_tools:build_defs.bzl", "gbbs_cc_library")

gbbs_cc_library(
    name = "scan",
    srcs = [
        "scan_helpers.cc",
    ],
    hdrs = [
        "scan.h",
        "scan_helpers.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//ligra",
        "//ligra:graph",
        "//ligra:macros",
        "//ligra:undirected_edge",
        "//ligra:vertex_subset",
        "//ligra/pbbslib:sparse_table",
        "//pbbslib:seq",
        "//pbbslib:stlalgs",
        "//pbbslib:utilities",
    ],
)
