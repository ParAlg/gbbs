licenses(["notice"])

cc_library(
    name = "SSWidestPath",
    hdrs = ["SSWidestPath.h"],
    deps = [
        "//gbbs",
        "//gbbs:bucket",
    ],
)

cc_binary(
    name = "SSWidestPath_main",
    srcs = ["SSWidestPath.cc"],
    deps = [":SSWidestPath"],
)

package(
    default_visibility = ["//visibility:public"],
)
