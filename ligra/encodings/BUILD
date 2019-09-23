cc_library(
  name = "byte",
  hdrs = ["byte.h"],
  deps = ["//ligra:macros"]
)

cc_library(
  name = "byte_pd",
  hdrs = ["byte_pd.h"],
  deps = [
  "//ligra:macros",
  "//pbbslib:binary_search",
  "//pbbslib:seq",
  "//pbbslib:sequence_ops",
  "//pbbslib:utilities",
  ]
)

cc_library(
  name = "byte_pd_amortized",
  hdrs = ["byte_pd_amortized.h"],
  deps = [
  "//ligra:bridge",
  ]
)


cc_library(
  name = "decoders",
  hdrs = ["decoders.h"],
  deps = [
  ":byte",
  ":byte_pd_amortized",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
