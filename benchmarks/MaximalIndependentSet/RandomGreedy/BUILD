cc_library(
  name = "MaximalIndependentSet",
  hdrs = ["MaximalIndependentSet.h"],
  deps = [
  "//ligra:ligra",
  "//ligra/pbbslib:sparse_table"
  ]
)

cc_binary(
  name = "MaximalIndependentSet_main",
  srcs = ["MaximalIndependentSet.cpp"],
  deps = [":MaximalIndependentSet"]
)

package(
  default_visibility = ["//visibility:public"],
)
