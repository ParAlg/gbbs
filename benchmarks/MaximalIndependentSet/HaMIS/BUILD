licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "HaMIS",
    hdrs = ["HaMIS.h"],
    deps = [
        "//gbbs",
    ],
)

cc_binary(
    name = "HaMIS_main",
    srcs = ["HaMIS.cc"],
    deps = [
        ":HaMIS",
        "//benchmarks/MaximalMatching/shared:decompress_to_csr",
        "//gbbs",
    ],
)
