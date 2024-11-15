load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
load("@bazel_tools//tools/build_defs/repo:git.bzl", "new_git_repository")
load("@bazel_tools//tools/cpp:cc_configure.bzl", "cc_configure")

cc_configure()

http_archive(
    name = "parlaylib",
    sha256 = "5ebff9645d6839c3c8fc6498b0089b69b39c9acec43556890596b51a420cc7f9",
    strip_prefix = "parlaylib-bazel-2024/include/",
    urls = ["https://github.com/ParAlg/parlaylib/archive/refs/tags/bazel-2024.tar.gz"],
)

http_archive(
    name = "googletest",
    sha256 = "b4870bf121ff7795ba20d20bcdd8627b8e088f2d1dab299a031c1034eddc93d5",
    strip_prefix = "googletest-release-1.11.0",
    urls = ["https://github.com/google/googletest/archive/release-1.11.0.tar.gz"],
)

#  path = external/parlaylib
#  url = https://github.com/ParAlg/parlaylib
#  branch = Bazel
