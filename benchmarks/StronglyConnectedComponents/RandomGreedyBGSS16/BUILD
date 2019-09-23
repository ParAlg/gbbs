cc_library(
  name = "StronglyConnectedComponents",
  hdrs = ["StronglyConnectedComponents.h"],
  deps = [
  "//ligra:ligra",
  "//ligra/pbbslib:resizable_table",
  "//pbbslib:random_shuffle"
  ]
)

cc_binary(
  name = "StronglyConnectedComponents_main",
  srcs = ["StronglyConnectedComponents.cpp"],
  deps = [":StronglyConnectedComponents"]
)

package(
  default_visibility = ["//visibility:public"],
)
