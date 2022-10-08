licenses(["notice"])

cc_library(
    name = "MaximalMatching",
    hdrs = ["MaximalMatching.h"],
    deps = [
        "//gbbs",
        "//gbbs/helpers:sparse_table",
        "//gbbs/helpers:speculative_for",
    ],
)

cc_binary(
    name = "MaximalMatching_main",
    srcs = ["MaximalMatching.cc"],
    deps = [
        ":MaximalMatching",
        "//gbbs",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
