cc_library(
  name = "wBFS",
  hdrs = ["wBFS.h"],
  deps = [
  "//ligra:ligra",
  "//ligra:bucket",
  ]
)

cc_binary(
  name = "wBFS_main",
  srcs = ["wBFS.cc"],
  deps = [":wBFS"]
)

package(
  default_visibility = ["//visibility:public"],
)
