licenses(["notice"])

cc_library(
    name = "byte",
    srcs = ["byte.cc"],
    hdrs = ["byte.h"],
    deps = ["//gbbs:macros"],
)

cc_library(
    name = "byte_pd",
    srcs = ["byte_pd.cc"],
    hdrs = ["byte_pd.h"],
    deps = [
        "//gbbs:bridge",
        "//gbbs:macros",
    ],
)

cc_library(
    name = "byte_pd_amortized",
    srcs = ["byte_pd_amortized.cc"],
    hdrs = ["byte_pd_amortized.h"],
    deps = [
        "//gbbs:bridge",
        "//gbbs:macros",
    ],
)

cc_library(
    name = "decoders",
    hdrs = ["decoders.h"],
    deps = [
        ":byte",
        ":byte_pd",
        ":byte_pd_amortized",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
