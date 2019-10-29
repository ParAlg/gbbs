cc_library(
  name = "Triangle",
  hdrs = ["Triangle.h"],
  deps = [
  "//ligra:ligra",
  "//pbbslib:sample_sort",
  "//pbbslib:monoid",
  ]
)

cc_binary(
  name = "Triangle_main",
  srcs = ["Triangle.cc"],
  deps = [":Triangle"]
)

package(
  default_visibility = ["//visibility:public"],
)
