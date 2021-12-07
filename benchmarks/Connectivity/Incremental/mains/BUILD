licenses(["notice"])

cc_binary(
    name = "jayanti_starting",
    srcs = ["jayanti.cc"],
    deps = [
        ":bench_utils",
        ":empty_starting_graph",
        ":subsample_starting_graph",
        ":uf_utils",
    ],
)

cc_binary(
    name = "jayanti_no_starting",
    srcs = ["jayanti.cc"],
    local_defines = ["EMPTY_STARTING_GRAPH"],
    deps = [
        ":bench_utils",
        ":empty_starting_graph",
        ":subsample_starting_graph",
        ":uf_utils",
    ],
)

cc_binary(
    name = "liutarjan_starting",
    srcs = ["liutarjan.cc"],
    deps = [
        ":bench_utils",
        ":empty_starting_graph",
        ":lt_utils",
        ":subsample_starting_graph",
    ],
)

cc_binary(
    name = "liutarjan_no_starting",
    srcs = ["liutarjan.cc"],
    local_defines = ["EMPTY_STARTING_GRAPH"],
    deps = [
        ":bench_utils",
        ":empty_starting_graph",
        ":lt_utils",
        ":subsample_starting_graph",
    ],
)

cc_binary(
    name = "shiloachvishkin_starting",
    srcs = ["shiloachvishkin.cc"],
    deps = [
        ":bench_utils",
        ":empty_starting_graph",
        ":subsample_starting_graph",
    ],
)

cc_binary(
    name = "shiloachvishkin_no_starting",
    srcs = ["shiloachvishkin.cc"],
    local_defines = ["EMPTY_STARTING_GRAPH"],
    deps = [
        ":bench_utils",
        ":empty_starting_graph",
        ":subsample_starting_graph",
    ],
)

cc_binary(
    name = "unite_starting",
    srcs = ["unite.cc"],
    deps = [
        ":bench_utils",
        ":empty_starting_graph",
        ":subsample_starting_graph",
        ":uf_utils",
    ],
)

cc_binary(
    name = "unite_no_starting",
    srcs = ["unite.cc"],
    local_defines = ["EMPTY_STARTING_GRAPH"],
    deps = [
        ":bench_utils",
        ":empty_starting_graph",
        ":subsample_starting_graph",
        ":uf_utils",
    ],
)

cc_binary(
    name = "unite_early_starting",
    srcs = ["unite_early.cc"],
    deps = [
        ":bench_utils",
        ":empty_starting_graph",
        ":subsample_starting_graph",
        ":uf_utils",
    ],
)

cc_binary(
    name = "unite_early_no_starting",
    srcs = ["unite_early.cc"],
    local_defines = ["EMPTY_STARTING_GRAPH"],
    deps = [
        ":bench_utils",
        ":empty_starting_graph",
        ":subsample_starting_graph",
        ":uf_utils",
    ],
)

cc_binary(
    name = "unite_nd_starting",
    srcs = ["unite_nd.cc"],
    deps = [
        ":bench_utils",
        ":empty_starting_graph",
        ":subsample_starting_graph",
        ":uf_utils",
    ],
)

cc_binary(
    name = "unite_nd_no_starting",
    srcs = ["unite_nd.cc"],
    local_defines = ["EMPTY_STARTING_GRAPH"],
    deps = [
        ":bench_utils",
        ":empty_starting_graph",
        ":subsample_starting_graph",
        ":uf_utils",
    ],
)

cc_binary(
    name = "unite_rem_cas_starting",
    srcs = ["unite_rem_cas.cc"],
    deps = [
        ":bench_utils",
        ":empty_starting_graph",
        ":subsample_starting_graph",
        ":uf_utils",
    ],
)

cc_binary(
    name = "unite_rem_cas_no_starting",
    srcs = ["unite_rem_cas.cc"],
    local_defines = ["EMPTY_STARTING_GRAPH"],
    deps = [
        ":bench_utils",
        ":empty_starting_graph",
        ":subsample_starting_graph",
        ":uf_utils",
    ],
)

cc_binary(
    name = "unite_rem_lock_starting",
    srcs = ["unite_rem_lock.cc"],
    deps = [
        ":bench_utils",
        ":empty_starting_graph",
        ":subsample_starting_graph",
        ":uf_utils",
    ],
)

cc_binary(
    name = "unite_rem_lock_no_starting",
    srcs = ["unite_rem_lock.cc"],
    local_defines = ["EMPTY_STARTING_GRAPH"],
    deps = [
        ":bench_utils",
        ":empty_starting_graph",
        ":subsample_starting_graph",
        ":uf_utils",
    ],
)

cc_library(
    name = "bench_utils",
    hdrs = ["bench_utils.h"],
    deps = [
        "//benchmarks/Connectivity/ConnectIt:framework",
        "//benchmarks/Connectivity/ConnectIt/mains:check",
        "//gbbs",
    ],
)

cc_library(
    name = "empty_starting_graph",
    hdrs = ["empty_starting_graph.h"],
    deps = [
    ],
)

cc_library(
    name = "lt_utils",
    hdrs = ["lt_utils.h"],
)

cc_library(
    name = "subsample_starting_graph",
    hdrs = ["subsample_starting_graph.h"],
    deps = [
        "//gbbs",
    ],
)

cc_library(
    name = "uf_utils",
    hdrs = ["uf_utils.h"],
    deps = [
        "//benchmarks/Connectivity:common",
    ],
)
