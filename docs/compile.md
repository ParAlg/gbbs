---
id: compile
title: Compile
---

## Compiling the Benchmarks

#### Using Bazel
You can build the benchmarks as follows:
```sh
# For Bazel:
$ bazel build //...  # compiles all benchmarks
```

Note that our bazel configuration file builds optimized binaries using
the`-c opt` mode by default. You can build a debug binary by
explicitly supplying the `-c dbg` option to bazel.

#### Using Make:
First set the appropriate environment variables. For example, if you
would like to compile with Cilk Plus first run `export CILK=1`. If no
options are set, the default is to use the homegrown scheduler.
After that, build using `make`.
```
$ cd benchmarks/BFS/NonDeterministicBFS  # go to a benchmark
$ make -j  # build the benchmark with all threads
```

## Cleaning the Build
The following commands cleans the directory:
```sh
# For Bazel:
$ bazel clean  # removes all executables

# For Make:
$ make clean  # removes executables for the current directory
```

## Other Build Options

#### Graphs with more than 2^32 Edges
To compile codes for graphs with more than 2^32 edges, the `LONG` command-line
parameter should be set. If the graph has more than 2^32 vertices, the
`EDGELONG` command-line parameter should be set. Note that the codes have not
been tested with more than 2^32 vertices, so if any issues arise
please mail [laxman@mit.edu](mailto:laxman@mit.edu).

#### Specifying Different Schedulers
To compile with the Cilk Plus scheduler instead of the Homegrown scheduler, use
the Bazel configuration `--config=cilk`. To compile using OpenMP instead, use
the Bazel configuration `--config=openmp`. To compile serially instead, use the
Bazel configuration `--config=serial`. (For the Makefiles, instead set the
environment variables `CILK`, `OPENMP`, or `SERIAL` respectively.)

