licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "wBFS",
    hdrs = ["wBFS.h"],
    deps = [
        "//gbbs",
        "//gbbs:bucket",
    ],
)

cc_binary(
    name = "wBFS_main",
    srcs = ["wBFS.cc"],
    deps = [":wBFS"],
)
