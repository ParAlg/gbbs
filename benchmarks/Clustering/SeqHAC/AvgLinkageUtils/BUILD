cc_library(
  name = "ClusteredGraph",
  hdrs = ["ClusteredGraph.h"],
  deps = [
  "//gbbs:gbbs",
  "@PAM//pam:pam",
  ]
)

cc_library(
  name = "HAC_configuration",
  hdrs = ["HAC_configuration.h"],
  deps = [
  "//gbbs:macros",
  ]
)

cc_library(
  name = "HeapBased",
  hdrs = ["HeapBased.h"],
  deps = [
  ":ClusteredGraph",
  "//gbbs:gbbs",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
