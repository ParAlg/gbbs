# GBBS: Graph Based Benchmark Suite

Organization
--------

This repository contains code for our SPAA paper "Theoretically Efficient
Parallel Graph Algorithms Can Be Fast and Scalable" (SPAA'18). It includes
implementations of the following parallel graph algorithms:

**Connectivity Problems**
* Low-Diameter Decomposition
* Connectivity
* Spanning Forest
* Biconnectivity
* Minimum Spanning Tree
* Strongly Connected Components

**Covering Problems**
* Coloring
* Maximal Matching
* Maximal Independent Set
* Approximate Set Cover

**Eigenvector Problems**
* PageRank

**Substructure Problems**
* Triangle Counting
* Approximate Densest Subgraph
* k-Core (coreness)

**Shortest Path Problems**
* Unweighted SSSP (Breadth-First Search)
* General Weight SSSP (Bellman-Ford)
* Integer Weight SSSP (Weighted Breadth-First Search)
* Single-Source Betweenness Centrality
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

Compiler:
* g++ &gt;= 5.3.0 with support for Cilk Plus
* g++ &gt;= 5.3.0 with pthread support (Homemade Scheduler)

Build system:
* [Bazel](https://docs.bazel.build/versions/master/install.html) 2.1.0

The default compilation uses a lightweight scheduler developed at CMU (Homemade)
for parallelism, which results in comparable performance to Cilk Plus. The
half-lengths for certain functions such as histogramming are lower using
Homemade, which results in better performance for codes like KCore.

Note: The Homemade scheduler was developed after our paper submission. For
reproducibility purposes, the codes should be compiled with Cilk Plus by adding
`--config=cilk`, although in our experience the times are usually faster using
Homemade.

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

To compile with the Cilk Plus scheduler instead of the Homegrown scheduler, use
the Bazel configuration `--config=cilk`. To compile using OpenMP instead, use
the Bazel configuration `--config=openmp`. To compile serially instead, use the
Bazel configuration `--config=serial`.

After setting the necessary environment variables:
```sh
$ bazel build --compilation_mode opt //...  # compiles all benchmarks
```

The following commands cleans the directory:
```sh
$ bazel clean  #removes all executables
```

Running code
-------
The applications take the input graph as input as well as an optional
flag "-s" to indicate a symmetric graph.  Symmetric graphs should be
called with the "-s" flag for better performance. For example:

```sh
$ bazel run --compilation_mode opt ///benchmarks/BFS/NonDeterministicBFS:BFS_main -- -s -src 10 ~/gbbs/inputs/rMatGraph_J_5_100
$ bazel run --compilation_mode opt ///benchmarks/IntegralWeightSSSP/JulienneDBS17:wBFS_main -- -s -w -src 15 ~/gbbs/inputs/rMatGraph_WJ_5_100
```

Note that the codes that compute single-source shortest paths (or centrality)
take an extra `-src` flag. The benchmark is run four times by default, and can
be changed by passing the `-rounds` flag followed by an integer indicating the
number of runs.

On NUMA machines, adding the command "numactl -i all " when running
the program may improve performance for large graphs. For example:

```sh
$ numactl -i all bazel run ...
```

Running code on compressed graphs
-----------

We make use of the bytePDA format in our benchmark, which is similar to the
parallelByte format of Ligra+, extended with additional functionality. We have
provided a converter utility which takes as input an uncompressed graph and
outputs a bytePDA graph. The converter can be used as follows:

```sh
bazel run --compilation_mode opt //utils:compressor -- -s -o ~/gbbs/inputs/rMatGraph_J_5_100.bytepda ~/gbbs/inputs/rMatGraph_J_5_100
bazel run --compilation_mode opt //utils:compressor -- -s -w -o ~/gbbs/inputs/rMatGraph_WJ_5_100.bytepda ~/gbbs/inputs/rMatGraph_WJ_5_100
```

After an uncompressed graph has been converted to the bytepda format,
applications can be run on it by passing in the usual command-line flags, with
an additional `-c` flag.

```sh
$ bazel run --compilation_mode opt //benchmarks/BFS/NonDeterministicBFS:BFS_main -- -s -c -src 10 ~/gbbs/inputs/rMatGraph_J_5_100.bytepda
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
