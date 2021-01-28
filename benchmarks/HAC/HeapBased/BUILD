cc_library(
  name = "HAC",
  hdrs = ["HAC.h"],
  deps = [
  "//gbbs:gbbs",
  ]
)

cc_binary(
  name = "HAC_main",
  srcs = ["HAC.cc"],
  deps = [":HAC"]
)

package(
  default_visibility = ["//visibility:public"],
)
