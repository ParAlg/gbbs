---
id: install
title: Install
---

## Requirements
#### System and Compiler Requirements
GBBS should work on Linux, Mac, and Windows machines, although primary
development and testing is done on Linux. The only requirement is a
relatively recent C++ compiler with C++17 support. Specifically, we
recommend:
* g++ >= 7.4.0 with support for Cilk Plus or
* g++ >= 7.4.0 with pthread support (Homemade Scheduler from
[ParlayLib](https://github.com/cmuparlay/parlaylib))

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

The default compilation uses a lightweight scheduler provided by
[ParlayLib](https://github.com/cmuparlay/parlaylib) for parallelism,
which results in comparable performance to Cilk Plus. The half-lengths
for certain functions such as histogramming are lower using Homemade,
which results in better performance for codes like KCore.


## Installation
Please download our repository from our
[github](https://www.github.com/ldhulipala/gbbs).

Note that if you are contributing to the repo and check out the repo
thorugh ssh, you can save yourself some nuisance re-typing your github
password every push by using
[ssh-agent + ssh-add](https://kb.iu.edu/d/aeww).
