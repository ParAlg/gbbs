licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "KCore",
    hdrs = ["KCore.h"],
    deps = [
        "//gbbs",
        "//gbbs:julienne",
    ],
)

cc_binary(
    name = "KCore_main",
    srcs = ["KCore.cc"],
    deps = [":KCore"],
)
