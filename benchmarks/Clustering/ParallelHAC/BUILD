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
    name = "HACConfiguration",
  hdrs = ["HACConfiguration.h"],
  deps = [
  "//gbbs:macros",
  ]
)

cc_binary(
  name = "ParallelHAC_main",
  srcs = ["ParallelHAC.cc"],
  deps = [
  ":ParallelHAC",
  ":HACConfiguration",
  ]
)


package(
  default_visibility = ["//visibility:public"],
)
