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
        "//benchmarks/SpanningForest:common",
        "//gbbs",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
