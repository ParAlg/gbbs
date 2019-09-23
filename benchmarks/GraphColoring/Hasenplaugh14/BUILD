cc_library(
  name = "Coloring",
  hdrs = ["Coloring.h"],
  deps = [
  "//ligra:ligra"
  "//pbbslib:random_shuffle"
  ]
)

cc_binary(
  name = "Coloring_main",
  srcs = ["Coloring.cpp"],
  deps = [":Coloring"]
)

package(
  default_visibility = ["//visibility:public"],
)
