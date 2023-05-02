licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "GraphColoring",
    hdrs = ["GraphColoring.h"],
    deps = [
        "//gbbs",
    ],
)

cc_binary(
    name = "GraphColoring_main",
    srcs = ["GraphColoring.cc"],
    deps = [":GraphColoring"],
)
