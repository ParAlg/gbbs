cc_library(
    name = "scan",
    srcs = [
        "scan.cc",
        "scan_helpers.cc",
    ],
    hdrs = [
        "scan.h",
        "scan_helpers.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":similarity_measure",
        ":utils",
        "//benchmarks/Connectivity/UnionFind:union_find_rules",
        "//gbbs:bridge",
        "//gbbs:graph",
        "//gbbs:macros",
        "//gbbs/pbbslib:sparse_table",
        "//pbbslib:binary_search",
        "//pbbslib:get_time",
        "//pbbslib:sample_sort",
        "//pbbslib:seq",
        "//pbbslib:utilities",
    ],
)

cc_library(
    name = "similarity_measure",
    srcs = [
        "similarity_measure.cc",
    ],
    hdrs = [
        "intersect.h",
        "similarity_measure.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//benchmarks/TriangleCounting/ShunTangwongsan15:Triangle",
        "//gbbs:bridge",
        "//gbbs:graph_mutation",
        "//gbbs:vertex",
        "//pbbslib:assert",
        "//pbbslib:binary_search",
        "//pbbslib:monoid",
        "//pbbslib:random",
        "//pbbslib:seq",
        "//pbbslib:sequence_ops",
        "//pbbslib:utilities",
    ],
)

cc_library(
    name = "utils",
    srcs = [
        "utils.cc",
    ],
    hdrs = [
        "utils.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//gbbs:bridge",
        "//gbbs:graph",
        "//gbbs:macros",
        "//pbbslib:integer_sort",
        "//pbbslib:monoid",
        "//pbbslib:seq",
    ],
)

cc_binary(
    name = "SCAN_main",
    srcs = ["run_scan.cc"],
    deps = [
        ":scan",
        ":similarity_measure",
        ":utils",
        "//gbbs",
    ],
)

cc_binary(
    name = "unweighted_experiments",
    srcs = ["unweighted_experiments.cc"],
    deps = [
        ":scan",
        ":similarity_measure",
        ":utils",
        "//gbbs",
        "//pbbslib:assert",
    ],
    linkopts = ["-lpython3.6m"],
)

cc_binary(
    name = "unweighted_experiments_serial",
    srcs = ["unweighted_experiments_serial.cc"],
    deps = [
        ":scan",
        ":similarity_measure",
        ":utils",
        "//gbbs",
        "//pbbslib:assert",
    ],
)

cc_binary(
    name = "weighted_experiments",
    srcs = ["weighted_experiments.cc"],
    deps = [
        ":scan",
        ":similarity_measure",
        ":utils",
        "//gbbs",
        "//pbbslib:assert",
    ],
    linkopts = ["-lpython3.6m"],
)

cc_binary(
    name = "weighted_experiments_serial",
    srcs = ["weighted_experiments_serial.cc"],
    deps = [
        ":scan",
        ":similarity_measure",
        ":utils",
        "//gbbs",
        "//pbbslib:assert",
    ],
)
