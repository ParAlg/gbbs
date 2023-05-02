licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "ApproximateSetCover",
    hdrs = ["ApproximateSetCover.h"],
    deps = [
        "//gbbs",
        "//gbbs:julienne",
    ],
)

cc_binary(
    name = "ApproximateSetCover_main",
    srcs = ["ApproximateSetCover.cc"],
    deps = [":ApproximateSetCover"],
)
