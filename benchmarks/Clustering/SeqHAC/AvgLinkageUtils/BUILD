licenses(["notice"])

cc_library(
    name = "ClusteredGraph",
    hdrs = ["ClusteredGraph.h"],
    deps = [
        "//gbbs",
        "@PAM//pam",
    ],
)

cc_library(
    name = "HeapBased",
    hdrs = ["HeapBased.h"],
    deps = [
        ":ClusteredGraph",
        "//gbbs",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
