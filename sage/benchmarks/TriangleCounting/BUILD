cc_library(
  name = "Triangle",
  hdrs = ["Triangle.h"],
  deps = [
  "//gbbs:gbbs",
  "//pbbslib:sample_sort",
  "//pbbslib:monoid",
  "//sage:sage",
  ]
)

cc_binary(
  name = "Triangle_main",
  srcs = ["Triangle.cc"],
  deps = [
    ":Triangle",
    "//sage:sage",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
