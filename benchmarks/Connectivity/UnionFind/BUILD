licenses(["notice"])

cc_library(
  name = "union_find_rules",
  hdrs = ["union_find_rules.h"],
  deps = [
  ":jayanti",
  "//gbbs:bridge",
  "//benchmarks/Connectivity:common",
  ]
)

cc_library(
  name = "jayanti",
  hdrs = ["jayanti.h"],
  deps = [
  "//gbbs:gbbs",
  "//benchmarks/Connectivity:common",
  ]
)

cc_library(
  name = "Connectivity",
  hdrs = ["Connectivity.h"],
  deps = [
  ":union_find_rules",
  "//gbbs:gbbs",
  "//gbbs/helpers:sparse_table",
  "//benchmarks/Connectivity:common",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
