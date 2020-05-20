cc_binary(
    name = "BarabasiAlbert",
    srcs = ["BarabasiAlbert.cc"],
    deps = [
        "BarabasiAlbert_library",
        "//gbbs",
        "//pbbslib/strings:string_basics",
    ],
)

cc_binary(
    name = "RMAT",
    srcs = ["RMAT.cc"],
    deps = [
        "RMAT_library",
        "//gbbs",
        "//pbbslib/strings:string_basics",
    ],
)

cc_library(
    name = "BarabasiAlbert_library",
    hdrs = ["barabasi_albert.h"],
    deps = [
        "//gbbs",
        "//gbbs:macros",
        "//pbbslib:random",
    ],
)

cc_library(
    name = "RMAT_library",
    hdrs = ["rmat.h"],
    deps = [
        "//gbbs",
        "//gbbs:macros",
    ],
)
