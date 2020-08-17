cc_library(
  name = "Triangle",
  hdrs = ["Triangle.h"],
  deps = [
  "//gbbs:gbbs",
  "//pbbslib:sample_sort",
  "//pbbslib:monoid",
  "//benchmarks/TriangleCounting/ShunTangwongsan15:Triangle",
  ]
)

cc_binary(
  name = "DynamicTriangle_main",
  srcs = ["Triangle.cc"],
  deps = [":Triangle"]
)

package(
  default_visibility = ["//visibility:public"],
)