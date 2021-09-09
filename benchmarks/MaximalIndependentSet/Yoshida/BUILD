licenses(["notice"])

cc_library(
  name = "reorder",
  hdrs = ["reorder.h"],
  deps = [
  "//gbbs:gbbs",
  ]
)

cc_library(
  name = "MaximalIndependentSet",
  hdrs = ["MaximalIndependentSet.h"],
  deps = [
  ":reorder",
  "//gbbs:gbbs",
  ]
)

cc_binary(
  name = "MaximalIndependentSet_main",
  srcs = ["MaximalIndependentSet.cc"],
  deps = [
  ":MaximalIndependentSet",
  "//gbbs:bridge",
  "//gbbs:gbbs",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
