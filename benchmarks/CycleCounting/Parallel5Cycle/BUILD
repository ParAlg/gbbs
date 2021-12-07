licenses(["notice"])

cc_library(
    name = "FiveCycle",
    hdrs = ["FiveCycle.h"],
    deps = [
        "//benchmarks/DegeneracyOrder/BarenboimElkin08:DegeneracyOrder",
        "//benchmarks/DegeneracyOrder/GoodrichPszona11:DegeneracyOrder",
        "//benchmarks/KCore/JulienneDBS17:KCore",
        "//gbbs",
        "//gbbs:graph",
        "//gbbs:graph_io",
    ],
)

cc_binary(
    name = "FiveCycle_main",
    srcs = ["FiveCycle.cc"],
    deps = [":FiveCycle"],
)

package(
    default_visibility = ["//visibility:public"],
)
