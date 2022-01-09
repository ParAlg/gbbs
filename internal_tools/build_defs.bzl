"""Custom Bazel build rules."""

def gbbs_cc_test(name, linkstatic = True, **kwargs):
    """Builds a C++ test. Analogous to the default `cc_test` Bazel rule."""
    # We set `linkstatic` to `True` by default (versus `cc_test` which sets it to
    # `False` by default). Rationale: the `-fvisibility=hidden` GCC flag makes
    # symbols private by default in shared libraries. With `linkstatic = False`,
    # then, unit tests will fail to compile when their tests involve non-public
    # symbols, which is most of the time. (In contrast, some symbols are public,
    # e.g., Python bindings that we purposely expose --- it would make sense to
    # set `linkstatic = False` for unit tests on these symbols.)

    native.cc_test(name = name, linkstatic = linkstatic, **kwargs)
