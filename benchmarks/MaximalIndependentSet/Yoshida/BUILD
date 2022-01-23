licenses(["notice"])

cc_library(
    name = "reorder",
    hdrs = ["reorder.h"],
    deps = [
        "//gbbs",
    ],
)

cc_library(
    name = "MaximalIndependentSet",
    hdrs = ["MaximalIndependentSet.h"],
    deps = [
        ":reorder",
        "//gbbs",
    ],
)

cc_binary(
    name = "MaximalIndependentSet_main",
    srcs = ["MaximalIndependentSet.cc"],
    deps = [
        ":MaximalIndependentSet",
        "//gbbs",
        "//gbbs:bridge",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
