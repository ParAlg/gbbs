cc_library(
    name = "KTruss",
    hdrs = [
        "KTruss.h",
        "truss_utils.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//gbbs:bridge",
        "//gbbs:bucket",
        "//gbbs:edge_map_reduce",
        "//gbbs:ligra",
        "//gbbs/pbbslib:dyn_arr",
        "//gbbs/pbbslib:sparse_table",
        "//pbbslib:sample_sort",
        "//pbbslib:utilities",
    ],
)

cc_binary(
    name = "KTruss_main",
    srcs = ["KTruss.cc"],
    deps = [":KTruss"],
)
