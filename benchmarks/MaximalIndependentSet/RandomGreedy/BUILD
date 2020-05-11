cc_library(
  name = "MaximalIndependentSet",
  hdrs = ["MaximalIndependentSet.h"],
  deps = [
  "//gbbs:ligra",
  "//gbbs/pbbslib:sparse_table"
  ]
)

cc_binary(
  name = "MaximalIndependentSet_main",
  srcs = ["MaximalIndependentSet.cc"],
  deps = [":MaximalIndependentSet"]
)

package(
  default_visibility = ["//visibility:public"],
)
