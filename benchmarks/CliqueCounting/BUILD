licenses(["notice"])

cc_library(
    name = "intersect",
    hdrs = ["intersect.h"],
    deps = [
        "//gbbs",
        "//gbbs:bucket",
        "//gbbs:edge_map_reduce",
    ],
)

cc_library(
    name = "induced_intersection",
    hdrs = ["induced_intersection.h"],
    deps = [
        ":intersect",
    ],
)

cc_library(
    name = "induced_neighborhood",
    hdrs = ["induced_neighborhood.h"],
    deps = [
        ":intersect",
    ],
)

cc_library(
    name = "induced_hybrid",
    hdrs = ["induced_hybrid.h"],
    deps = [
        ":intersect",
    ],
)

cc_library(
    name = "induced_split",
    hdrs = ["induced_split.h"],
    deps = [
        ":induced_hybrid",
        ":intersect",
    ],
)

cc_library(
    name = "peel",
    hdrs = ["peel.h"],
    deps = [
        ":induced_hybrid",
        ":intersect",
        "//gbbs",
        "//gbbs:bucket",
        "//gbbs:edge_map_reduce",
        "//gbbs/helpers:sparse_table",
    ],
)

cc_library(
    name = "relabel",
    hdrs = ["relabel.h"],
    deps = [
        "//gbbs",
        "//gbbs:graph",
    ],
)

cc_library(
    name = "Clique",
    hdrs = ["Clique.h"],
    deps = [
        ":induced_hybrid",
        ":induced_intersection",
        ":induced_neighborhood",
        ":induced_split",
        ":intersect",
        ":peel",
        ":relabel",
        "//benchmarks/DegeneracyOrder/BarenboimElkin08:DegeneracyOrder",
        "//benchmarks/DegeneracyOrder/GoodrichPszona11:DegeneracyOrder",
        "//benchmarks/KCore/JulienneDBS17:KCore",
        "//benchmarks/TriangleCounting/ShunTangwongsan15:Triangle",
        "//gbbs",
        "//gbbs:bucket",
        "//gbbs:edge_map_reduce",
        "//gbbs/helpers:assert",
    ],
)

cc_binary(
    name = "Clique_main",
    srcs = ["Clique.cc"],
    deps = [
        ":Clique",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
