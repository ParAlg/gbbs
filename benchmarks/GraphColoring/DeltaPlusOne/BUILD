cc_library(
  name = "DynamicDeltaPlusOneColoring",
  hdrs = ["DynamicDeltaPlusOneColoring.h"],
  deps = [
  "//gbbs:gbbs",
  "//pbbslib:random",
  ]
)

cc_binary(
  name = "DynamicDeltaPlusOneColoring_main",
  srcs = ["DynamicDeltaPlusOneColoring.cc"],
  deps = [":DynamicDeltaPlusOneColoring"]
)

package(
  default_visibility = ["//visibility:public"],
)
