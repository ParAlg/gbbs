licenses(["notice"])

cc_library(
    name = "union_find_rules",
    hdrs = ["union_find_rules.h"],
    deps = [
        ":jayanti",
        "//benchmarks/Connectivity:common",
        "//gbbs:bridge",
    ],
)

cc_library(
    name = "jayanti",
    hdrs = ["jayanti.h"],
    deps = [
        "//benchmarks/Connectivity:common",
        "//gbbs",
    ],
)

cc_library(
    name = "Connectivity",
    hdrs = ["Connectivity.h"],
    deps = [
        ":union_find_rules",
        "//benchmarks/Connectivity:common",
        "//gbbs",
        "//gbbs/helpers:sparse_table",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
