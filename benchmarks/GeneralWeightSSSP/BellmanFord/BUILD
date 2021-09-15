licenses(["notice"])

cc_library(
    name = "BellmanFord",
    hdrs = ["BellmanFord.h"],
    deps = [
        "//gbbs",
    ],
)

cc_binary(
    name = "BellmanFord_main",
    srcs = ["BellmanFord.cc"],
    deps = [":BellmanFord"],
)

package(
    default_visibility = ["//visibility:public"],
)
