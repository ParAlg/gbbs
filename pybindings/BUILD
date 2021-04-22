package(default_visibility = ["//visibility:public"])

load("@pybind11_bazel//:build_defs.bzl", "pybind_extension", "pybind_library")

pybind_extension(
  name = "gbbs_lib",
  srcs = ["gbbs_lib.cc"],
  deps = [
      "//gbbs:compressed_vertex",
      "//gbbs:graph",
      "//gbbs:graph_io",
      "//gbbs:gbbs",
      "//gbbs:vertex",
      "//gbbs:vertex_subset",
      "//benchmarks/BFS/NonDeterministicBFS:BFS",
      "//benchmarks/Biconnectivity/TarjanVishkin:Biconnectivity",
      "//benchmarks/Clustering/SeqHAC:HAC_api",
      "//benchmarks/Connectivity/WorkEfficientSDB14:Connectivity",
      "//benchmarks/GeneralWeightSSSP/BellmanFord:BellmanFord",
      "//benchmarks/MinimumSpanningForest/Boruvka:MinimumSpanningForest",
      "//benchmarks/StronglyConnectedComponents/RandomGreedyBGSS16:StronglyConnectedComponents",
      "//benchmarks/CoSimRank:CoSimRank",
      "//benchmarks/PageRank:PageRank",
      "//benchmarks/KCore/JulienneDBS17:KCore",
  ],
)