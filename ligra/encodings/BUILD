load("//internal_tools:build_defs.bzl", "gbbs_cc_library")

gbbs_cc_library(
  name = "byte",
  hdrs = ["byte.h"],
  srcs = ["byte.cc"],
  deps = ["//ligra:macros"]
)

gbbs_cc_library(
  name = "byte_pd",
  hdrs = ["byte_pd.h"],
  deps = [
  "//ligra:macros",
  "//pbbslib:binary_search",
  "//pbbslib:seq",
  "//pbbslib:sequence_ops",
  "//pbbslib:utilities",
  ]
)

gbbs_cc_library(
  name = "byte_pd_amortized",
  hdrs = ["byte_pd_amortized.h"],
  srcs = ["byte_pd_amortized.cc"],
  deps = [
  "//ligra:bridge",
  "//ligra:macros",
  ]
)


gbbs_cc_library(
  name = "decoders",
  hdrs = ["decoders.h"],
  deps = [
  ":byte",
  ":byte_pd_amortized",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
