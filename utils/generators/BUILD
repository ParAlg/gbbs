load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_binary(
    name = "BarabasiAlbert",
    srcs = ["BarabasiAlbert.cc"],
    deps = [
        "BarabasiAlbert_library",
        "//ligra",
        "//ligra:bridge",
        "//pbbslib/strings:string_basics",
    ],
)

gbbs_cc_binary(
    name = "RMAT",
    srcs = ["RMAT.cc"],
    deps = [
        "RMAT_library",
        "//ligra",
        "//ligra:bridge",
        "//pbbslib/strings:string_basics",
    ],
)

gbbs_cc_library(
    name = "BarabasiAlbert_library",
    hdrs = ["barabasi_albert.h"],
    deps = [
        "//ligra",
        "//ligra:macros",
        "//pbbslib:random",
    ],
)

gbbs_cc_library(
    name = "RMAT_library",
    hdrs = ["rmat.h"],
    deps = [
        "//ligra",
        "//ligra:macros",
    ],
)
