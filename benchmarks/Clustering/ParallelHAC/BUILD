cc_library(
  name = "ClusteredGraph",
  hdrs = ["ClusteredGraph.h"],
  deps = [
  "//gbbs:gbbs",
  "@PAM//pam:pam",
  ]
)

cc_library(
  name = "ParallelHAC",
  hdrs = ["ParallelHAC.h"],
  deps = [
  ":ClusteredGraph",
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

cc_binary(
  name = "ParallelHAC_main",
  srcs = ["ParallelHAC.cc"],
  deps = [
  ":ParallelHAC",
  ":HAC_configuration",
  ]
)


package(
  default_visibility = ["//visibility:public"],
)
