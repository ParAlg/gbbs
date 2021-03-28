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
Library) and place the installation so that the `setvars.sh` script is located
at `/home/ubuntu/intel/oneapi/setvars.sh` (or otherwise modify the SCAN makefile accordingly).

The Bazel build file does not work with MKL, so build with `make` instead.
