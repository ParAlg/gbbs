"""Custom Bazel build rules."""

# C++ compiler flags compatible with all supported compilers.
CXX_FLAGS = [
    "-std=c++17",
    "-march=native",
    "-mcx16",  # enable 16-byte compare-and-swap
    "-fvisibility=hidden",
    "-fvisibility-inlines-hidden",
    "-DLONG",  # use 8 byte vertex identifiers
    "-DAMORTIZEDPD",  # use amortized_bytepd encoding scheme for compressed graphs
    "-DUSEMALLOC",
    "-Wall",
    "-Wextra",
    "-Wcast-qual",
    "-Wno-unused-parameter",
    "-Wpointer-arith",
    "-Wvla",
    # TODO(tomtseng): clean up code to enable -Wshadow (and remove -Wshadow=local
    # from GCC_CXX_FLAGS) and -Wmissing-declarations
]

def _get_compiler_specific_cxx_flags():
    """Get compiler flags for a specific C++ compiler (GCC, Clang)."""
    CLANG_CXX_FLAGS = []
    GCC_CXX_FLAGS = [
        "-Wshadow=local",
    ]
    return select(
        {
            # Assume that macOS users use Clang and Linux users use GCC.
            "//internal_tools:apple": CLANG_CXX_FLAGS,
            "//internal_tools:linux": GCC_CXX_FLAGS,
            "//conditions:default": [],
        }
    )


def _get_parallelism_cxx_flags():
    """Returns C++ compiler flags specific to the parallelism framework."""
    SERIAL_CXX_FLAGS = []
    CILK_CXX_FLAGS = [
        "-DCILK",
        "-fcilkplus",
    ]
    OPENMP_CXX_FLAGS = [
        "-DOPENMP",
        "-fopenmp",
    ]
    HOMEGROWN_SCHEDULER_CXX_FLAGS = [
        "-DHOMEGROWN",
        "-pthread",
    ]
    return select(
        {
            "//internal_tools:serial": SERIAL_CXX_FLAGS,
            "//internal_tools:cilk": CILK_CXX_FLAGS,
            "//internal_tools:openmp": OPENMP_CXX_FLAGS,
            # Use homegrown scheduler by default.
            "//conditions:default": HOMEGROWN_SCHEDULER_CXX_FLAGS,
        }
    )


def _get_parallelism_linker_flags():
    """Returns C++ linker flags specific to the parallelism framework."""
    SERIAL_LINKER_FLAGS = []
    CILK_LINKER_FLAGS = ["-lcilkrts"]
    OPENMP_LINKER_FLAGS = ["-fopenmp"]
    HOMEGROWN_SCHEDULER_LINKER_FLAGS = []
    return select(
        {
            "//internal_tools:serial": SERIAL_LINKER_FLAGS,
            "//internal_tools:cilk": CILK_LINKER_FLAGS,
            "//internal_tools:openmp": OPENMP_LINKER_FLAGS,
            # Use homegrown scheduler by default.
            "//conditions:default": HOMEGROWN_SCHEDULER_LINKER_FLAGS,
        }
    )


def _get_copts():
    """Returns all C++ compiler flags."""
    return CXX_FLAGS + _get_compiler_specific_cxx_flags() + _get_parallelism_cxx_flags()


def _get_linkopts():
    """Returns all C++ linker flags."""
    return _get_parallelism_linker_flags()


def gbbs_cc_binary(name, copts=[], linkopts=[], **kwargs):
    """Builds a C++ binary. Analogous to the default `cc_binary` Bazel rule."""

    native.cc_binary(
        name=name,
        copts=_get_copts() + copts,
        linkopts=_get_linkopts() + linkopts,
        **kwargs
    )


def gbbs_cc_library(name, copts=[], linkopts=[], **kwargs):
    """Builds a C++ library. Analogous to the default `cc_library` Bazel
    rule."""

    native.cc_library(
        name=name,
        copts=_get_copts() + copts,
        linkopts=_get_linkopts() + linkopts,
        **kwargs
    )


def gbbs_cc_test(name, copts=[], linkopts=[], linkstatic=True, **kwargs):
    """Builds a C++ test. Analogous to the default `cc_test` Bazel rule."""
    # We set `linkstatic` to `True` by default (versus `cc_test` which sets it to
    # `False` by default). Rationale: the `-fvisibility=hidden` GCC flag makes
    # symbols private by default in shared libraries. With `linkstatic = False`,
    # then, unit tests will fail to compile when their tests involve non-public
    # symbols, which is most of the time. (In contrast, some symbols are public,
    # e.g., Python bindings that we purposely expose --- it would make sense to
    # set `linkstatic = False` for unit tests on these symbols.)

    native.cc_test(
        name=name,
        copts=_get_copts() + copts,
        linkopts=_get_linkopts() + linkopts,
        linkstatic=linkstatic,
        **kwargs
    )
