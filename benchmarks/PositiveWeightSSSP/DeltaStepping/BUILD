licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "DeltaStepping",
    hdrs = ["DeltaStepping.h"],
    deps = [
        "//gbbs",
        "//gbbs:bucket",
    ],
)

cc_binary(
    name = "DeltaStepping_main",
    srcs = ["DeltaStepping.cc"],
    deps = [":DeltaStepping"],
)
