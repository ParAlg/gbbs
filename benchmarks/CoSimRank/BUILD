licenses(["notice"])

cc_library(
    name = "CoSimRank",
    hdrs = ["CoSimRank.h"],
    deps = [
        "//gbbs",
    ],
)

cc_binary(
    name = "CoSimRank_main",
    srcs = ["CoSimRank.cc"],
    deps = [":CoSimRank"],
)

package(
    default_visibility = ["//visibility:public"],
)
