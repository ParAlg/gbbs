# Index-based SCAN

## About

This implementation of SCAN is a parallel version of the index-based SCAN
algorithm introduced in "Efficient structural graph clustering: an index-based
approach" by Wen et al.

To invoke the implementation, see the `Index` class in `scan.h`.

## Additional notes

Define the `SCAN_DETAILED_TIMES` macro in order to output more detailed timings.

## MKL

To use matrix multiplication on dense graphs, install Intel MKL (Math Kernel
Library). The Bazel build file is not set up to work with MKL, so build with `make`
instead. Source MKL's `setvars.sh` file before running `make`.
