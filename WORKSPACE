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


# Creates a repository rule for the system python headers.
# pybind11.BUILD depends on this repository rule to detect your python configuration
load("//third_party/pybind11_bazel:python_configure.bzl", "python_configure")
python_configure(name = "local_config_python")

# Create pybind11 external repository
# If using another pybind11 version:
# Use tar URL of desired version, change strip_prefix to your version "pybind11-x.x.x",
# Supply correct sha256 for your version.
http_archive(
    name = "pybind11",
    build_file = "@//third_party/pybind11_bazel:pybind11.BUILD",
    sha256 = "1eed57bc6863190e35637290f97a20c81cfe4d9090ac0a24f3bbf08f265eb71d",
    strip_prefix = "pybind11-2.4.3",
    url = "https://github.com/pybind/pybind11/archive/v2.4.3.tar.gz",
)

# @rules_python repository, used to create python build targets
http_archive(
    name = "rules_python",
    sha256 = "aa96a691d3a8177f3215b14b0edc9641787abaaa30363a080165d06ab65e1161",
    url = "https://github.com/bazelbuild/rules_python/releases/download/0.0.1/rules_python-0.0.1.tar.gz",
)

# Currently does nothing, futureproofs your core Python rule dependencies.
load("@rules_python//python:repositories.bzl", "py_repositories")
py_repositories()

# Pulls in dependencies needed to use the python packaging rules.
load("@rules_python//python:pip.bzl", "pip_repositories")
pip_repositories()
