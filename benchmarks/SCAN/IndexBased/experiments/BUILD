licenses(["notice"])

#cc_binary(
#    name = "run_gbbs_experiments",
#    srcs = ["run_gbbs_experiments.cc"],
#    linkopts = ["-lpython3.6m"],
#    deps = [
#        "//benchmarks/SCAN/IndexBased:scan",
#        "//benchmarks/SCAN/IndexBased:similarity_measure",
#        "//benchmarks/SCAN/IndexBased:utils",
#        "//gbbs",
#        "//helpers:assert",
#    ],
#)

cc_binary(
    name = "remove_singleton_vertices",
    srcs = ["remove_singleton_vertices.cc"],
    deps = [
        "//gbbs/helpers:parse_command_line",
    ],
)
