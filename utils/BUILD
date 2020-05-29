load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_binary(
    name = "compressor",
    srcs = ["compressor.cc"],
    deps = [
        "//ligra",
        "//ligra:io",
        "//ligra:parse_command_line",
        "//pbbslib:random",
        "//pbbslib:utilities",
    ],
)

gbbs_cc_binary(
    name = "converter",
    srcs = ["converter.cc"],
    deps = [
        ":converter_library",
    ],
)

gbbs_cc_binary(
    name = "random_reorder",
    srcs = ["random_reorder.cc"],
    deps = [
        "//ligra",
        "//ligra:bridge",
        "//pbbslib:random_shuffle",
        "//pbbslib/strings:string_basics",
    ],
)

gbbs_cc_binary(
    name = "to_edge_list",
    srcs = ["to_edge_list.cc"],
    deps = [
        "//ligra",
        "//ligra:bridge",
        "//pbbslib:random_shuffle",
        "//pbbslib/strings:string_basics",
    ],
)

gbbs_cc_binary(
    name = "snap_converter",
    srcs = ["snap_converter.cc"],
    deps = [
        "//ligra:graph_io",
        "//ligra:parse_command_line",
    ],
)

gbbs_cc_library(
    name = "converter_library",
    hdrs = [
        "converter.h",
        "to_char_arr.h",
    ],
    deps = ["//ligra"],
)
