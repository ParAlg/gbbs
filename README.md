# GBBS: Graph Based Benchmark Suite

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
* Mininmum Spanning Tree
* Strongly Connected Components
* Approximate Set Cover
* Triangle Counting
* Weighted Breadth-First Search

The code for these applications is located in the `benchmark` directory. The
implementations are based on the Ligra/Ligra+/Julienne graph processing
frameworks. The framework code is located in the `src` directory.

The codes used here are still in development, and we plan to add more
applications/benchmarks. We currently include the following extra codes:

* Densest Subgraph (the (2+\epsilon)-approximation from Bahmani et al.)
* KTruss

If you use our work, please cite our paper:

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

The default compilation uses CILK Plus. We also support a lightweight scheduler
developed at CMU (Homemade), which results in comparable performance to CILK.
The half-lengths for certain functions such as histogramming are lower using
Homemade, which results in better performance for codes like KCore.

Note: The Homemade scheduler was developed after our paper submission. For
reproducibility purposes, the codes should be compiled with CILK Plus, although
in our experience the times are usually faster using Homemade.

The benchmark supports both uncompressed and compressed graphs. The uncompressed
format is identical to the uncompressed format in Ligra. The compressed format
is the parallelByte format used in Ligra+, with some additional functionality to
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


Running code
-------
The applications take the input graph as input as well as an optional
flag "-s" to indicate a symmetric graph.  Symmetric graphs should be
called with the "-s" flag for better performance. For example:

```
$ ./BFS -s -src 10 ../inputs/rMatGraph_J_5_100
$ ./BellmanFord -s -src 15 ../inputs/rMatGraph_WJ_5_100
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


Ongoing Work
--------
We are currently working on:

* writing unit-tests

* porting over utilities for converting large graphs to the bytepda format, and
  adding random edge weights to the graphs.
