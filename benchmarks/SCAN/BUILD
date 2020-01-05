cc_library(
    name = "scan",
    srcs = ["symmetric_2d_table.cc"],
    hdrs = [
        "scan.h",
        "symmetric_2d_table.h",
    ],
    deps = [
        "//ligra:macros",
        "//utils:assert",
    ],
)

cc_binary(
    name = "run_scan",
    srcs = ["run_scan.cc"],
    deps = [
        ":scan",
        "//ligra",
        "//utils:assert",
    ],
)
