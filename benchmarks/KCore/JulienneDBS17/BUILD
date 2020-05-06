cc_library(
  name = "KCore",
  hdrs = ["KCore.h"],
  deps = [
  "//ligra:gbbs",
  "//ligra:julienne",
  ]
)

cc_binary(
  name = "KCore_main",
  srcs = ["KCore.cc"],
  deps = [":KCore"]
)

package(
  default_visibility = ["//visibility:public"],
)
