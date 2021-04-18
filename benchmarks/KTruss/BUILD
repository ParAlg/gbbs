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
        "//gbbs:gbbs",
        "//gbbs/pbbslib:dyn_arr",
        "//gbbs/pbbslib:sparse_table",
    ],
)

cc_binary(
    name = "KTruss_main",
    srcs = ["KTruss.cc"],
    deps = [":KTruss"],
)
