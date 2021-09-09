licenses(["notice"])

cc_library(
  name = "reorder",
  hdrs = ["reorder.h"],
  deps = [
    "//gbbs:gbbs",
  ]
)

cc_library(
  name = "yoshida_matching_utils",
  hdrs = ["yoshida_matching_utils.h"],
  deps = []
)

cc_library(
  name = "MaximalMatching",
  hdrs = ["MaximalMatching.h"],
  deps = [
  ":reorder",
  ":yoshida_matching_utils",
  "//gbbs:bridge",
  "//gbbs:gbbs",
  ]
)

cc_binary(
  name = "MaximalMatching_main",
  srcs = ["MaximalMatching.cc"],
  deps = [
  ":MaximalMatching",
  "//gbbs:bridge",
  "//gbbs:gbbs",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
