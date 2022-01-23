licenses(["notice"])

cc_binary(
    name = "compressor",
    srcs = ["compressor.cc"],
    deps = [
        "//gbbs",
        "//gbbs:io",
        "//gbbs/encodings:byte_pd_amortized",
        "//gbbs/helpers:parse_command_line",
    ],
)

cc_binary(
    name = "converter",
    srcs = ["converter.cc"],
    deps = [
        ":converter_library",
    ],
)

cc_binary(
    name = "random_reorder",
    srcs = ["random_reorder.cc"],
    deps = [
        "//gbbs",
    ],
)

cc_binary(
    name = "to_edge_list",
    srcs = ["to_edge_list.cc"],
    deps = [
        "//gbbs",
    ],
)

cc_binary(
    name = "snap_converter",
    srcs = ["snap_converter.cc"],
    deps = [
        "//gbbs:graph_io",
        "//gbbs/helpers:parse_command_line",
    ],
)

cc_library(
    name = "converter_library",
    hdrs = [
        "converter.h",
        "to_char_arr.h",
    ],
    deps = ["//gbbs"],
)
