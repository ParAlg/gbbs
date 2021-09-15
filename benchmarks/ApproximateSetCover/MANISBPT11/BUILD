licenses(["notice"])

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

package(
    default_visibility = ["//visibility:public"],
)
