cc_library(
  name = "bench_utils",
  hdrs = ["bench_utils.h"],
  deps = []
)

cc_binary(
  name = "EdgeIncremental",
  srcs = ["EdgeIncremental.cc"],
  deps = [":bench_utils",
          "//benchmarks/Connectivity/WorkEfficientSDB14:Connectivity",
          "//benchmarks/Connectivity/UnionFind:Connectivity",
          "//benchmarks/Connectivity/LabelPropagation:Connectivity",
          "//benchmarks/Connectivity/BFSCC:Connectivity",
          "//benchmarks/Connectivity/ShiloachVishkin:Connectivity",
          "//benchmarks/Connectivity/Framework:framework",
          "//benchmarks/Connectivity:common",
          "//ligra:graph_mutation",
          "//pbbslib/strings:string_basics"
         ]
)

package(
  default_visibility = ["//visibility:public"],
)
