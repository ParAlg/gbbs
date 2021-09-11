load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
load("@bazel_tools//tools/build_defs/repo:git.bzl", "new_git_repository")

load("@bazel_tools//tools/cpp:cc_configure.bzl", "cc_configure")
cc_configure()

local_repository(
    name = "parlaylib",
    path = "external/parlaylib/include",
)

local_repository(
    name = "PAM",
    path = "external/PAM/include",
)

http_archive(
  name = "googletest",
  urls = ["https://github.com/google/googletest/archive/release-1.11.0.tar.gz"],
  strip_prefix = "googletest-release-1.11.0",
  sha256 = "b4870bf121ff7795ba20d20bcdd8627b8e088f2d1dab299a031c1034eddc93d5",
)

# pybind bazel bindings
PYBIND11_BAZEL_COMMIT = "656fd69e83e80ccef8a3e3935d1ff4f604332e81"
http_archive(
  name = "pybind11_bazel",
  strip_prefix = "pybind11_bazel-%s" % PYBIND11_BAZEL_COMMIT,
  sha256 = "b17dbbb648001996d06cce7e058174687c933bb16d87b74b340565d5311bb286",
  urls = ["https://github.com/ldhulipala/pybind11_bazel/archive/%s.zip" % PYBIND11_BAZEL_COMMIT],
)
# pybind.
PYBIND11_COMMIT = "fe755dce12766820a99eefbde32d6ceb0a828ca8"
http_archive(
  name = "pybind11",
  build_file = "@pybind11_bazel//:pybind11.BUILD",
  strip_prefix = "pybind11-%s" % PYBIND11_COMMIT,
  sha256 = "5702350060e965043609a5115576c598ef3262a7367f063aebfeaa7961f2cbfd",
  urls = ["https://github.com/pybind/pybind11/archive/%s.tar.gz" % PYBIND11_COMMIT],
)

load("@pybind11_bazel//:python_configure.bzl", "python_configure")
python_configure(name = "local_config_python")
