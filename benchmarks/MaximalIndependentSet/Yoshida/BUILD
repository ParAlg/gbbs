load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "reorder",
  hdrs = ["reorder.h"],
  deps = [
  "//ligra:bridge",
  "//ligra:ligra",
  ]
)

gbbs_cc_library(
  name = "MaximalIndependentSet",
  hdrs = ["MaximalIndependentSet.h"],
  deps = [
  ":reorder",
  "//ligra:bridge",
  "//ligra:ligra",
  ]
)

gbbs_cc_binary(
  name = "MaximalIndependentSet_main",
  srcs = ["MaximalIndependentSet.cc"],
  deps = [
  ":MaximalIndependentSet",
  "//ligra:bridge",
  "//ligra:ligra",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
