licenses(["notice"])

cc_library(
    name = "PushRelabel",
    hdrs = ["PushRelabel.h"],
    deps = ["//gbbs"],
)

cc_binary(
    name = "PushRelabel_main",
    srcs = ["PushRelabel.cc"],
    deps = [":PushRelabel"],
)

package(
    default_visibility = ["//visibility:public"],
)
