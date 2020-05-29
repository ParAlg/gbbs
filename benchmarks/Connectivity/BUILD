load("//internal_tools:build_defs.bzl", "gbbs_cc_library")

gbbs_cc_library(
  name = "common",
  hdrs = [
  "common.h",
  "connectit.h",
  ],
  srcs = [
  "common.cc",
  "connectit.cc",
  ],
  deps = [
  "//ligra/pbbslib:atomic_counter",
  "//pbbslib:assert",
  ],
)

package(
  default_visibility = ["//visibility:public"],
)
