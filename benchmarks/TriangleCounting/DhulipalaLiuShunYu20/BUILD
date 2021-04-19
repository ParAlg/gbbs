cc_library(
  name = "sparse_table",
  hdrs = ["sparse_table.h"],
  deps = [
    "//gbbs:bridge"
  ]
)

cc_library(
  name = "tomb_table",
  hdrs = ["tomb_table.h"],
  deps = [
    "//gbbs:bridge"
  ]
)

cc_library(
  name = "shared",
  hdrs = ["shared.h"],
  deps = [
    "//gbbs:gbbs",
  ]
)

cc_library(
  name = "dynamic_graph",
  hdrs = ["dynamic_graph.h"],
  deps = [
    "//gbbs:gbbs",
    ":shared",
    ":sparse_table",
    ":tomb_table",
  ]
)

cc_library(
  name = "preprocess",
  hdrs = ["preprocess.h"],
  deps = [
    "//gbbs:gbbs",
    ":dynamic_graph",
  ]
)

cc_library(
  name = "rebalancing",
  hdrs = ["rebalancing.h"],
  deps = [
    "//gbbs:gbbs",
    ":shared",
    ":sparse_table",
    ":tomb_table",
    ":dynamic_graph",
    ":preprocess",
    "//benchmarks/TriangleCounting/ShunTangwongsan15:Triangle",
  ]
)

cc_library(
  name = "benchmark",
  hdrs = ["benchmark.h"],
  deps = [
    "//gbbs:gbbs",
    ":shared",
    ":rebalancing",
  ]
)

cc_library(
  name = "Triangle",
  hdrs = ["Triangle.h"],
  deps = [
  "//gbbs:gbbs",
  "//benchmarks/TriangleCounting/ShunTangwongsan15:Triangle",
  ":shared",
  ":benchmark",
  ":preprocess",
  ":dynamic_graph",
  ":rebalancing",
  ]
)


cc_library(
  name = "makkar_graph",
  hdrs = ["makkar_graph.h"],
  deps = [
    "//gbbs:gbbs",
  ]
)

cc_library(
  name = "Makkar",
  hdrs = ["makkar.h"],
  deps = [
    "//gbbs:gbbs",
    ":makkar_graph",
  ]
)

cc_binary(
  name = "DynamicTriangle_main",
  srcs = ["Triangle.cc"],
  deps = [":Triangle", ":Makkar"]
)

package(
  default_visibility = ["//visibility:public"],
)
