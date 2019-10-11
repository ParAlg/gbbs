load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")
load("@bazel_tools//tools/build_defs/repo:git.bzl", "new_git_repository")

GRAPHSETINTER_BUILD = """
genrule(
  name = "intersection_genrule",
  srcs = glob(["src/set_operation.cpp", "src/set_operation.hpp"]),
  outs = ["set_operation.o"],
  cmd = "cd external/graphsetinter/src; make set_operation.o; cd ../../..; cp external/graphsetinter/src/set_operation.o $(location set_operation.o)",
)

cc_library(
  name = "intersectiondeps",
  hdrs = glob(["src/*.hpp"]),
  visibility = ["//visibility:public"],
)

cc_library(
  name = "intersection",
  srcs = ["set_operation.o"],
  hdrs = ["src/set_operation.hpp"],
  deps = ["//:intersectiondeps"],
  visibility = ["//visibility:public"],
)
"""

new_git_repository(
  name = "graphsetinter",
  remote = "https://github.com/jeshi96/GraphSetIntersection.git",
  commit = "27339f4cfec1f989d39dff67bd39bee26748de7b",
  build_file_content = GRAPHSETINTER_BUILD,
)

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
