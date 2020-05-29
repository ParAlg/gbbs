load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "MaximalIndependentSet",
  hdrs = ["MaximalIndependentSet.h"],
  deps = [
  "//ligra:ligra",
  "//ligra/pbbslib:sparse_table"
  ]
)

gbbs_cc_binary(
  name = "MaximalIndependentSet_main",
  srcs = ["MaximalIndependentSet.cc"],
  deps = [":MaximalIndependentSet"]
)

package(
  default_visibility = ["//visibility:public"],
)
