licenses(["notice"])

cc_library(
  name = "FiveCycle",
  hdrs = ["FiveCycle.h"],
  deps = [
  "//gbbs:gbbs",
  "//gbbs:graph",
  "//benchmarks/DegeneracyOrder/BarenboimElkin08:DegeneracyOrder",
  "//benchmarks/DegeneracyOrder/GoodrichPszona11:DegeneracyOrder",
  "//benchmarks/KCore/JulienneDBS17:KCore"
  ]
)

cc_binary(
  name = "FiveCycle_main",
  srcs = ["FiveCycle.cc"],
  deps = [":FiveCycle"]
)

package(
  default_visibility = ["//visibility:public"],
)
