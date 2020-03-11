cc_binary(
    name = "BarabasiAlbert",
    srcs = ["BarabasiAlbert.cc"],
    deps = [
        "BarabasiAlbert_library",
        "//ligra",
        "//ligra:bridge",
        "//pbbslib/strings:string_basics",
    ],
)

cc_binary(
    name = "RMAT",
    srcs = ["RMAT.cc"],
    deps = [
        "RMAT_library",
        "//ligra",
        "//ligra:bridge",
        "//pbbslib/strings:string_basics",
    ],
)

cc_library(
    name = "BarabasiAlbert_library",
    hdrs = ["barabasi_albert.h"],
    deps = [
        "//ligra",
        "//ligra:macros",
        "//pbbslib:random",
    ],
)

cc_library(
    name = "RMAT_library",
    hdrs = ["rmat.h"],
    deps = [
        "//ligra",
        "//ligra:macros",
    ],
)
