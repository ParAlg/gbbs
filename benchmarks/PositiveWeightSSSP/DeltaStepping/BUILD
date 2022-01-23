licenses(["notice"])

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

package(
    default_visibility = ["//visibility:public"],
)
