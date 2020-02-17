cc_library(
  name = "reorder",
  hdrs = ["reorder.h"],
  deps = [
  "//ligra:bridge",
  "//ligra:ligra",
  ]
)

cc_library(
  name = "MaximalIndependentSet",
  hdrs = ["MaximalIndependentSet.h"],
  deps = [
  ":reorder",
  "//ligra:bridge",
  "//ligra:ligra",
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
