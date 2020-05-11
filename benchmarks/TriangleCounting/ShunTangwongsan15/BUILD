cc_library(
  name = "Triangle",
  hdrs = ["Triangle.h"],
  deps = [
  "//benchmarks/DegeneracyOrder/BarenboimElkin08:DegeneracyOrder",
  "//benchmarks/DegeneracyOrder/GoodrichPszona11:DegeneracyOrder",
  "//benchmarks/KCore/JulienneDBS17:KCore",
  "//gbbs:ligra",
  "//pbbslib:sample_sort",
  "//pbbslib:monoid",
  ]
)

cc_binary(
  name = "Triangle_main",
  srcs = ["Triangle.cc"],
  deps = [":Triangle"]
)

package(
  default_visibility = ["//visibility:public"],
)
