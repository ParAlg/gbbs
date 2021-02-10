# Parallel Batch-Dynamic k-Core Decomposition and Low Out-Degree Orientation

Organization
--------

This repository contains the implementations for our
approximate k-core algorithms.

`benchmarks/EdgeOrientation/` contains the implementations for our batch-dynamic
approximate k-core algorithm. Specifically, `ParallelLDS` (PLDS) contains the
implementation for our parallel batch-dynamic
![equation](https://latex.codecogs.com/gif.latex?%282&plus;%5Cdelta%29)-approximate
k-core decomposition algorithm, while `LDS` contains the baseline sequential
implementation of the same algorithm.

`benchmarks/KCore/ApproximateKCore` contains the implementation for our static
![equation](https://latex.codecogs.com/gif.latex?%282&plus;%5Cdelta%29)-approximate
k-core decomposition algorithm.

Compilation
--------

Compiler:
* g++ &gt;= 7.4.0 with support for Cilk Plus
* g++ &gt;= 7.4.0 with pthread support (Parlay Scheduler)

Build system:
* Make

To build:
```sh
$ cd benchmarks/EdgeOrientation/ParallelLDS  # go to a benchmark
$ make
```

The following commands cleans the directory:
```sh
$ make clean  # removes executables for the current directory
```

Most optionality from the [Graph Based Benchmark Suite (GBBS)](https://github.com/ParAlg/gbbs) and
[ParlayLib](https://github.com/cmuparlay/parlaylib) apply. In particular, to compile benchmarks for graphs with
more than 2^32 edges, the `LONG` command-line parameter should be set.

### Using Standard Library Allocator

The default compilation uses a scalable memory allocator developed at CMU
(Parlay). However, it is also possible to use the C++ standard library,
by adding the compile definition `-DPARLAY_USE_STD_ALLOC`.

### Using Cilk, OpenMP, or TBB

The default compilation uses a lightweight scheduler developed at CMU (Parlay)
for parallelism, which results in comparable performance to Cilk Plus.
However, it is also possible to use a Cilk, OpenMP, or Thread Building
Blocks with our implementations instead of the Parlay scheduler; simply add
the appropriate compile definition as below.

```
-DPARLAY_CILK
-DPARLAY_OPENMP
-DPARLAY_TBB
```

### Graph Format

The applications take dynamic graphs as input (in all three benchmarks) or 
static graphs as input (in the static benchmark).

For the static graph format, we support the adjacency graph format used by
[Graph Based Benchmark Suite (GBBS)](https://github.com/ParAlg/gbbs).

For the dynamic graph format, the file should contain one dynamic edge per line
(in order), with a `+` or a `-` indicating an insertion or a deletion
respectively, followed by the two vertex IDs. We assume all
vertices are in the range [0, n]. For instance, a single line of the file
could be `+ 0 1`, to insert edge (0, 1). The file should be represented
as plain text.