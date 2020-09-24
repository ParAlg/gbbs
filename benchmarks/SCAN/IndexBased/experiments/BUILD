cc_binary(
    name = "unweighted_experiments",
    srcs = ["unweighted_experiments.cc"],
    linkopts = ["-lpython3.6m"],
    deps = [
        "//benchmarks/SCAN/IndexBased:scan",
        "//benchmarks/SCAN/IndexBased:similarity_measure",
        "//benchmarks/SCAN/IndexBased:utils",
        "//gbbs",
        "//pbbslib:assert",
    ],
)

cc_binary(
    name = "weighted_experiments",
    srcs = ["weighted_experiments.cc"],
    linkopts = ["-lpython3.6m"],
    deps = [
        "//benchmarks/SCAN/IndexBased:scan",
        "//benchmarks/SCAN/IndexBased:similarity_measure",
        "//benchmarks/SCAN/IndexBased:utils",
        "//gbbs",
        "//pbbslib:assert",
    ],
)
