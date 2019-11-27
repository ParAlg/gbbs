cc_library(
  name = "framework",
  hdrs = ["framework.h"],
  deps = [":sampling",
          "//benchmarks/Connectivity/WorkEfficientSDB14:Connectivity",
          "//benchmarks/Connectivity/UnionFind:Connectivity",
          "//benchmarks/Connectivity/LabelPropagation:Connectivity",
          "//benchmarks/Connectivity/ShiloachVishkin:Connectivity",
          "//benchmarks/Connectivity/LiuTarjan:Connectivity",
          "//benchmarks/Connectivity/BFSCC:Connectivity",
          "//benchmarks/Connectivity:common"]
)

cc_library(
  name = "sampling",
  hdrs = ["sampling.h"],
  deps = [":utils"]
)


cc_library(
  name = "utils",
  hdrs = ["utils.h"],
  deps = []
)

cc_binary(
  name = "ConnectivityBenchmark",
  srcs = ["ConnectivityBenchmark.cc"],
  deps = [":framework",
          "//benchmarks/Connectivity/WorkEfficientSDB14:Connectivity",
          "//benchmarks/Connectivity/BFSCC:Connectivity",
          "//ligra:graph_mutation",
          "//pbbslib/strings:string_basics"
         ]
)

package(
  default_visibility = ["//visibility:public"],
)
