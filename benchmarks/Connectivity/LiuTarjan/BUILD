cc_library(
  name = "liu_tarjan_rules",
  hdrs = ["liu_tarjan_rules.h"],
  deps = [
  "//benchmarks/Connectivity:common",
  ]
)

cc_library(
  name = "Connectivity",
  hdrs = ["Connectivity.h"],
  deps = [
  ":liu_tarjan_rules",
  "//gbbs:gbbs",
  "//pbbslib:random",
  "//benchmarks/Connectivity:common",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
