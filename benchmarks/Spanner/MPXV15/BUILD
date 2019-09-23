cc_library(
  name = "Spanner",
  hdrs = ["Spanner.h"],
  deps = [
  "//ligra:ligra",
  "//ligra/pbbslib:sparse_table"
  "//ligra/pbbslib:dyn_arr"
  ]
)

cc_binary(
  name = "Spanner_main",
  srcs = ["Spanner.cpp"],
  deps = [":Spanner"]
)

package(
  default_visibility = ["//visibility:public"],
)
