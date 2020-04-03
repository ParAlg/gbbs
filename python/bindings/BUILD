package(default_visibility = ["//visibility:public"])

load("@pybind11_bazel//:build_defs.bzl", "pybind_extension", "pybind_library")

pybind_extension(
  name = "gbbs_lib",
  srcs = ["gbbs_lib.cc"],
  deps = [
      "//ligra:compressed_vertex",
      "//ligra:graph",
      "//ligra:graph_io",
      "//ligra:ligra",
      "//ligra:vertex",
      "//ligra:vertex_subset",
  ],
)
