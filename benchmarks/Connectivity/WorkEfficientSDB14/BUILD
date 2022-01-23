licenses(["notice"])

cc_library(
    name = "Connectivity",
    hdrs = ["Connectivity.h"],
    deps = [
        "//benchmarks/Connectivity:common",
        "//benchmarks/LowDiameterDecomposition/MPX13:LowDiameterDecomposition",
        "//gbbs",
        "//gbbs:contract",
    ],
)

cc_binary(
    name = "Connectivity_main",
    srcs = ["Connectivity.cc"],
    deps = [":Connectivity"],
)

package(
    default_visibility = ["//visibility:public"],
)
