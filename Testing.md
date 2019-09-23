# Running Tests and Creating New Tests


### Adding new tests

You can add a new c++ test using googletest (example).


To hook the test up using Bazel, add the following to the BUILD-file, adding in
any extra dependencies that might be required for the test.

```
cc_test(
  name = "example_test",
  srcs = ["example_test.cpp"],
  copts = ["-Iexternal/gtest/include"],
  deps = [
  "@googletest//:gtest_main",
  ]
)
```


