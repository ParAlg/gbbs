licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "graph_filter",
    hdrs = ["graph_filter.h"],
    deps = [
        ":bitset_managers",
        "//gbbs:flags",
        "//gbbs:graph",
    ],
)

cc_library(
    name = "utils",
    hdrs = ["utils.h"],
    deps = [
        "//gbbs:macros",
    ],
)

cc_library(
    name = "bitset_managers",
    hdrs = ["bitset_managers.h"],
    deps = [
        ":bitset",
        ":utils",
        "//gbbs:flags",
        "//gbbs:macros",
        "//gbbs/encodings:byte_pd_amortized",
    ],
)

cc_library(
    name = "bitset",
    hdrs = ["bitset.h"],
    deps = [
        "//gbbs:macros",
    ],
)
