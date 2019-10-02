# GBBS: Graph Based Benchmark Suite [![Build Status](https://travis-ci.com/ldhulipala/gbbs.svg?branch=master)](https://travis-ci.com/ldhulipala/gbbs)

Organization
--------

This repository contains code for our SPAA paper "Theoretically Efficient
Parallel Graph Algorithms Can Be Fast and Scalable" (SPAA'18). It includes
implementations of the following parallel graph algorithms:

* Single-Source Betweenness Centrality
* Bellman-Ford
* Breadth-First Search
* Biconnectivity
* Connectivity
* Coloring
* KCore
* Low-Diameter Decomposition
* Maximal Matching
* Maximal Independent Set
* Minimum Spanning Tree
* Strongly Connected Components
* Approximate Set Cover
* Triangle Counting
* Weighted Breadth-First Search
* Approximate Densest Subgraph
* Spanning Forest
* PageRank
* Single-Source Widest Path
* k-Spanner 

The code for these applications is located in the `benchmark` directory. The
implementations are based on the Ligra/Ligra+/Julienne graph processing
frameworks. The framework code is located in the `src` directory.

The codes used here are still in development, and we plan to add more
applications/benchmarks. We currently include the following extra codes,
which are part of ongoing work.

* experimental/KTruss

If you use our work, please cite our [paper](https://arxiv.org/abs/1805.05208):

```
@inproceedings{dhulipala2018theoretically,
  author    = {Laxman Dhulipala and
               Guy E. Blelloch and
               Julian Shun},
  title     = {Theoretically Efficient Parallel Graph Algorithms Can Be Fast and
               Scalable},
  booktitle = {ACM Symposium on Parallelism in Algorithms and Architectures (SPAA)},
  year      = {2018},
}
```

Compilation
--------

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
format is identical to the uncompressed format in Ligra. The compressed format,
called bytepd_amortized (bytepda) is similar to the parallelByte format used in
Ligra+, with some additional functionality to support efficiently packs,
filters, and other operations over neighbor lists.

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

Running code
-------
The applications take the input graph as input as well as an optional
flag "-s" to indicate a symmetric graph.  Symmetric graphs should be
called with the "-s" flag for better performance. For example:

```
$ ./BFS -s -src 10 ../inputs/rMatGraph_J_5_100
$ ./wBFS -s -w -src 15 ../inputs/rMatGraph_WJ_5_100
```

Note that the codes that compute single-source shortest paths (or centrality)
take an extra `-src` flag. The benchmark is run four times by default, and can
be changed by passing the `-rounds` flag followed by an integer indicating the
number of runs.

On NUMA machines, adding the command "numactl -i all " when running
the program may improve performance for large graphs. For example:

```
$ numactl -i all ./BFS -s <input file>
```

Running code on compressed graphs
-----------

We make use of the bytePDA format in our benchmark, which is similar to the
parallelByte format of Ligra+, extended with additional functionality. We have
provided a converter utility which takes as input an uncompressed graph and
outputs a bytePDA graph. The converter can be used as follows:

```
./compressor -s -o ../inputs/rMatGraph_J_5_100.bytepda ../inputs/rMatGraph_J_5_100
./compressor -s -w -o ../inputs/rMatGraph_WJ_5_100.bytepda ../inputs/rMatGraph_WJ_5_100
```

After an uncompressed graph has been converted to the bytepda format,
applications can be run on it by passing in the usual command-line flags, with
an additional `-c` flag.

```
$ ./BFS -s -c -src 10 ../inputs/rMatGraph_J_5_100.bytepda
$ ./wBFS -s -w -c -src 15 ../inputs/rMatGraph_WJ_5_100.bytepda
```

When processing large compressed graphs, using the `-m` command-line flag can
help if the file is already in the page cache, since the compressed graph data
can be mmap'd. Application performance will be affected if the file is not
already in the page-cache. We have found that using `-m` when the compressed
graph is backed by SSD results in a slow first-run, followed by fast subsequent
runs.


Input Formats
-----------
We support the adjacency graph format used by the [Problem Based Benchmark
suite](http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html)
and [Ligra](https://github.com/jshun/ligra).

The adjacency graph format starts with a sequence of offsets one for each
vertex, followed by a sequence of directed edges ordered by their source vertex.
The offset for a vertex i refers to the location of the start of a contiguous
block of out edges for vertex i in the sequence of edges. The block continues
until the offset of the next vertex, or the end if i is the last vertex. All
vertices and offsets are 0 based and represented in decimal. The specific format
is as follows:

```
AdjacencyGraph
<n>
<m>
<o0>
<o1>
...
<o(n-1)>
<e0>
<e1>
...
<e(m-1)>
```

This file is represented as plain text.

Weighted graphsare represented in the weighted adjacnecy graph format. The file
should start with the string "WeightedAdjacencyGraph". The m edges weights
should be stored after all of the edge targets in the .adj file.
