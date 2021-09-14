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
  "//benchmarks/SpanningForest:common",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
