cc_library(
    name = "KTruss",
    hdrs = [
        "KTruss.h",
        "truss_utils.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//ligra",
        "//ligra:bridge",
        "//ligra:bucket",
        "//ligra:edge_map_reduce",
        "//ligra/pbbslib:dyn_arr",
        "//ligra/pbbslib:sparse_table",
        "//pbbslib:sample_sort",
        "//pbbslib:utilities",
    ],
)

cc_binary(
    name = "KTruss_main",
    srcs = ["KTruss.cc"],
    deps = [":KTruss"],
)
