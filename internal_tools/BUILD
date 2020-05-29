# Platform configuration labels.

config_setting(
    name = "apple",
    constraint_values = ["@bazel_tools//platforms:osx"],
)

config_setting(
    name = "linux",
    constraint_values = ["@bazel_tools//platforms:linux"],
)

# Parallelism configuration labels.

config_setting(
    name = "serial",
    define_values = {"parallel": "no"},
)

config_setting(
    name = "cilk",
    define_values = {"parallel": "cilk"},
)

config_setting(
    name = "openmp",
    define_values = {"parallel": "openmp"},
)
