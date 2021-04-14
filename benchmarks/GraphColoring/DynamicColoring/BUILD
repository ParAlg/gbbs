cc_library(
  name = "DynamicColoring",
  hdrs = ["DynamicColoring.h"],
  deps = [
  "//gbbs:gbbs",
  "//pbbslib:random",
  ]
)

cc_binary(
  name = "DynamicColoring_main",
  srcs = ["DynamicColoring.cc"],
  deps = [":DynamicColoring"]
)

package(
  default_visibility = ["//visibility:public"],
)
