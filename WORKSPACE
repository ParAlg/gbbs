load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")
load("@bazel_tools//tools/build_defs/repo:git.bzl", "new_git_repository")

SIMDINTER_BUILD = """
genrule(
  name = "intersection_genrule",
  srcs = glob(["src/intersection.cpp", "include/intersection.h"]),
  outs = ["intersection.o"],
  cmd = "cd external/simdinter; make intersection.o; cd ../..; cp external/simdinter/intersection.o $(location intersection.o)",
)

cc_library(
  name = "intersectiondeps",
  hdrs = glob(["include/*.h"]),
  visibility = ["//visibility:public"],
)

cc_library(
  name = "intersection",
  srcs = ["intersection.o"],
  hdrs = ["include/intersection.h"],
  deps = ["//:intersectiondeps"],
  visibility = ["//visibility:public"],
)
"""

new_git_repository(
  name = "simdinter",
  remote = "https://github.com/lemire/SIMDCompressionAndIntersection.git",
  commit = "f002db1d47f252dd17daa8206e3ebbbeee9e4d9b",
  build_file_content = SIMDINTER_BUILD,
)

git_repository(
  name = "googletest",
  remote = "https://github.com/google/googletest",
  tag = "release-1.8.1",
)
