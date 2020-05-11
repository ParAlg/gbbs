package(
  default_visibility = ["//visibility:public"],
)

cc_test(
    name = "graph_io_test",
    srcs = ["graph_io_test.cc"],
    deps = [
        "//gbbs:graph_io",
        "//gbbs:graph_test_utils",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "graph_test",
    srcs = ["graph_test.cc"],
    deps = [
        "//gbbs:graph",
        "//pbbslib:seq",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "undirected_edge_test",
    srcs = ["undirected_edge_test.cc"],
    deps = [
        "//gbbs:undirected_edge",
        "@googletest//:gtest_main",
    ],
)
