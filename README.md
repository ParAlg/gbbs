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


