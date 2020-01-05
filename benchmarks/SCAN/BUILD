cc_library(
    name = "scan",
    srcs = ["scan.cc"],
    hdrs = ["scan.h"],
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
