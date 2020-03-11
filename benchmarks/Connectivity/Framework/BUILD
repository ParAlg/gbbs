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
  deps = []
)

package(
  default_visibility = ["//visibility:public"],
)
