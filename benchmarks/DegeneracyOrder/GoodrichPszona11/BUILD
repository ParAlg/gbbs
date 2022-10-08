licenses(["notice"])

cc_library(
    name = "DegeneracyOrder",
    hdrs = ["DegeneracyOrder.h"],
    deps = [
        "//gbbs",
        "//gbbs:bucket",
        "//gbbs:edge_map_reduce",
    ],
)

cc_binary(
    name = "DegeneracyOrder_main",
    srcs = ["DegeneracyOrder.cc"],
    deps = [":DegeneracyOrder"],
)

package(
    default_visibility = ["//visibility:public"],
)
