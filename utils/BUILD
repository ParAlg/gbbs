cc_binary(
    name = "compressor",
    srcs = ["compressor.cc"],
    deps = [
        "//gbbs",
        "//gbbs:io",
        "//gbbs:parse_command_line",
        "//gbbs/encodings:byte_pd_amortized",
        "//pbbslib:random",
        "//pbbslib:utilities",
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
        "//pbbslib:random_shuffle",
        "//pbbslib/strings:string_basics",
    ],
)

cc_binary(
    name = "to_edge_list",
    srcs = ["to_edge_list.cc"],
    deps = [
        "//gbbs",
        "//pbbslib:random_shuffle",
        "//pbbslib/strings:string_basics",
    ],
)

cc_binary(
    name = "snap_converter",
    srcs = ["snap_converter.cc"],
    deps = [
        "//gbbs:graph_io",
        "//gbbs:parse_command_line",
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
