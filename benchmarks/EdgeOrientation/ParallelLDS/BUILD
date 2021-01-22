cc_library(
  name = "LDS",
  hdrs = ["LDS.h"],
  deps = [
  "//gbbs:gbbs",
  ]
)

cc_binary(
  name = "LDS_main",
  srcs = ["LDS.cc"],
  deps = [":LDS"]
)

package(
  default_visibility = ["//visibility:public"],
)
