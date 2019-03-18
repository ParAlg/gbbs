---
title: "Requirements"
layout: archive
permalink: /benchmarks/requirements/
comments: false
mathjax: true
sidebar:
  nav: "benchmarks"
---

* g++ &gt;= 5.3.0 with support for Cilk Plus
* g++ &gt;= 5.3.0 with pthread support (Homemade Scheduler)

The default compilation uses Cilk Plus. We also support a lightweight scheduler
developed at CMU (Homemade), which results in comparable performance to Cilk.
The half-lengths for certain functions such as histogramming are lower using
Homemade, which results in better performance for codes like KCore.

Note: The Homemade scheduler was developed after our paper submission. For
reproducibility purposes, the codes should be compiled with Cilk Plus, although
in our experience the times are usually faster using Homemade.

The benchmark supports both uncompressed and compressed graphs. The uncompressed
[format](/benchmarks/formats/) is identical to the uncompressed format in Ligra.
The compressed format, called bytepd_amortized (bytepda) is similar to the
parallelByte format used in Ligra+, with some additional functionality to
support efficiently packs, filters, and other operations over neighbor lists.

To compile codes for graphs with more than 2^32 edges, the `LONG` command-line
parameter should be set. If the graph has more than 2^32 vertices, the
`EDGELONG` command-line parameter should be set. Note that the codes have not
been tested with more than 2^32 vertices, so if any issues arise please contact
[Laxman Dhulipala](mailto:ldhulipa@cs.cmu.edu).

To compile using the Homemade scheduler the `HOMEMADE` command-line parameter
should be set. If it is unset, the Cilk Plus scheduler is used by default.

After setting the necessary environment variables:
```
$ make -j  #compiles the benchmark with all threads
```

The following commands cleans the directory:
```
$ make clean  #removes all executables
```
