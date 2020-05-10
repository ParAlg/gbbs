cc_library(
  name = "reorder",
  hdrs = ["reorder.h"],
  deps = [
  "//ligra:gbbs",
  ]
)

cc_library(
  name = "MaximalIndependentSet",
  hdrs = ["MaximalIndependentSet.h"],
  deps = [
  ":reorder",
  "//ligra:gbbs",
  ]
)

cc_binary(
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
