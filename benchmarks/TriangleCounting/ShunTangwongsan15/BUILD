load("//internal_tools:build_defs.bzl", "gbbs_cc_binary", "gbbs_cc_library")

gbbs_cc_library(
  name = "Triangle",
  hdrs = ["Triangle.h"],
  deps = [
  "//benchmarks/DegeneracyOrder/BarenboimElkin08:DegeneracyOrder",
  "//benchmarks/DegeneracyOrder/GoodrichPszona11:DegeneracyOrder",
  "//benchmarks/KCore/JulienneDBS17:KCore",
  "//ligra:ligra",
  "//pbbslib:sample_sort",
  "//pbbslib:monoid",
  ]
)

gbbs_cc_binary(
  name = "Triangle_main",
  srcs = ["Triangle.cc"],
  deps = [":Triangle"]
)

package(
  default_visibility = ["//visibility:public"],
)
