load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
    name = "KTruss",
    hdrs = [
        "KTruss.h",
        "truss_utils.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//ligra:bridge",
        "//ligra:bucket",
        "//ligra:edge_map_reduce",
        "//ligra:ligra",
        "//ligra/pbbslib:dyn_arr",
        "//ligra/pbbslib:sparse_table",
        "//pbbslib:sample_sort",
        "//pbbslib:utilities",
    ],
)

gbbs_cc_binary(
    name = "KTruss_main",
    srcs = ["KTruss.cc"],
    deps = [":KTruss"],
)
