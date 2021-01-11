cc_library(
  name = "byte",
  hdrs = ["byte.h"],
  srcs = ["byte.cc"],
  deps = ["//gbbs:macros"]
)

cc_library(
  name = "byte_pd",
  hdrs = ["byte_pd.h"],
  srcs = ["byte_pd.cc"],
  deps = [
  "//gbbs:bridge",
  "//gbbs:macros",
  ]
)

cc_library(
  name = "byte_pd_amortized",
  hdrs = ["byte_pd_amortized.h"],
  srcs = ["byte_pd_amortized.cc"],
  deps = [
  "//gbbs:bridge",
  "//gbbs:macros",
  ]
)


cc_library(
  name = "decoders",
  hdrs = ["decoders.h"],
  deps = [
  ":byte",
  ":byte_pd",
  ":byte_pd_amortized",
  ]
)

package(
  default_visibility = ["//visibility:public"],
)
