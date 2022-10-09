licenses(["notice"])

cc_library(
    name = "KCoreOrdering",
    hdrs = ["KCoreOrdering.h"],
    deps = [
        "//gbbs",
        "//gbbs:julienne",
        "//benchmarks/KCore/JulienneDBS17:KCore",
        "//benchmarks/Reordering:reorder_utils",
    ],
)

cc_binary(
    name = "KCoreOrdering_main",
    srcs = ["KCoreOrdering.cc"],
    deps = [":KCoreOrdering"],
)

package(
    default_visibility = ["//visibility:public"],
)
