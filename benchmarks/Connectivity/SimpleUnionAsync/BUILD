licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "Connectivity",
    hdrs = ["Connectivity.h"],
    deps = [
        "//benchmarks/Connectivity:common",
        "//gbbs",
    ],
)

cc_binary(
    name = "Connectivity_main",
    srcs = ["Connectivity.cc"],
    deps = [":Connectivity"],
)
