licenses(["notice"])

cc_library(
    name = "Biconnectivity",
    hdrs = ["Biconnectivity.h"],
    deps = [
        "//benchmarks/Connectivity/SimpleUnionAsync:Connectivity",
        "//gbbs",
    ],
)

cc_binary(
    name = "Biconnectivity_main",
    srcs = ["Biconnectivity.cc"],
    deps = [
        ":Biconnectivity",
        "//gbbs/helpers:sparse_additive_map",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
