cc_library(
  name = "fullyDynamic",
  hdrs = ["Hypergraph.h", 
    "GraphScheduler.h", 
    "HypergraphCoreDecomp.h"],
)

cc_library(
  name = "KCore",
  hdrs = ["FullyDynamic.h"],
  deps = [
  "//gbbs:gbbs",
  ":fullyDynamic"
  ]
)

cc_binary(
  name = "KCore_main",
  srcs = ["KCore.cc"],
  deps = [":KCore"]
)

package(
  default_visibility = ["//visibility:public"],
)
