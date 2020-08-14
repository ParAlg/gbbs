cc_library(
  name = "Triangle",
  hdrs = ["Triangle.h"],
  deps = [
  "//gbbs:gbbs",
  "//pbbslib:sample_sort",
  "//pbbslib:monoid",
  "//benchmarks/TriangleCounting/ShunTangwongsan15:Triangle"
  ":sparse_table",
  ":dynamic_graph",
  ":benchmark",
  ":preprocess",
  ":rebalancing",
  ":shared",
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