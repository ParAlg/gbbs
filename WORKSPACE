load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
load("@bazel_tools//tools/build_defs/repo:git.bzl", "new_git_repository")

load("@bazel_tools//tools/cpp:cc_configure.bzl", "cc_configure")
cc_configure()

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
  commit = "ca7e88ae14b804f1ec6889b989544053d842b1c9",
  shallow_since = "1571170990 -0400",
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
  shallow_since = "1562290323 -0400",
  build_file_content = SIMDINTER_BUILD,
)

http_archive(
  name = "googletest",
  urls = ["https://github.com/google/googletest/archive/release-1.10.0.tar.gz"],
  strip_prefix = "googletest-release-1.10.0",
  sha256 = "9dc9157a9a1551ec7a7e43daea9a694a0bb5fb8bec81235d8a1e6ef64c716dcb",
)

load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")

# Add support for Boost using https://github.com/nelhage/rules_boost {{{
git_repository(
    name = "com_github_nelhage_rules_boost",
    commit = "9f9fb8b2f0213989247c9d5c0e814a8451d18d7f",
    remote = "https://github.com/nelhage/rules_boost",
    shallow_since = "1570056263 -0700",
)

load("@com_github_nelhage_rules_boost//:boost/boost.bzl", "boost_deps")
boost_deps()
# }}} end adding support for Boost
