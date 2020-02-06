cc_library(
  name = "FiveCycle",
  hdrs = ["FiveCycle.h"],
  deps = [
  "//ligra:ligra",
  "//pbbslib:random_shuffle"
  "//pbbslib:sample_sort",
  ]
)

cc_binary(
  name = "FiveCycle_main",
  srcs = ["FiveCycle.cc"],
  deps = [":FiveCycle"]
)

package(
  default_visibility = ["//visibility:public"],
)
