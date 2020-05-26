---
id: install
title: Install and Compile
---

## Requirements
#### System and Compiler Requirements
GBBS should work on Linux, Mac, and Windows machines, although primary
development and testing is done on Linux. The only requirement is a
relatively recent C++ compiler with C++17 support. Specifically, we
recommend:
* g++ >= 7.4.0 with pthread support (using our custom scheduler)

If you would like to use Cilk for scheduling, you need:
* g++ >= 7.4.0 with support for Cilk Plus

You should also be able to use GBBS in conjunction with a version of
LLVM supporting C++17 that has support for Cilk linguistic features,
which can be found [here](http://cilk.mit.edu/).

#### Build System Requirements
We support two methods for building our code:
* [Bazel](https://docs.bazel.build/versions/master/install.html). This is the primary build method that we use.
* Make. 

Note that Make is a secondary method for building and tinkering with
GBBS that we support for users who do not wish to install Bazel.
However, all of our performance numbers are collected through binaries
built by Bazel, and so Bazel should be used if using the benchmarks in
a performance-critical setting.


## Installation
Please download our repository from our
[github](https://www.github.com/ldhulipala/gbbs).

Note that if you are contributing to the repo and check out the repo
thorugh ssh, you can save yourself some nuisance re-typing your github
password every push by using
[ssh-agent + ssh-add](https://kb.iu.edu/d/aeww).

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
First set the appropriate environment variables, e.g., first run
`export CILK=1` to compile with Cilk Plus.
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

