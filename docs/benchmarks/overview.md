---
id: overview
title: Overview
description: Overview of GBBS Benchmarks
permalink: /docs/benchmarks/overview
---

This page specifies the graph problems considered in GBBS. Accessing
one of the benchmark pages below provides the following information:
* Formal input-output specification of the problem.
* Algorithms solving the problem. Specifically, for each algorithm we provide:
  * Pointers to the code on github implementing the algorithm
  * Work and depth bounds of the algorithm
  * References to the papers/sources for the algorithm
* Examples of how to compile and run an implementation on a graph input.

We note that we are continuing to work on new benchmarks, and new
implementations of existing problems, and so the data here may only be
partially complete (please check our [github
repo](https://www.github.com/ParAlg/gbbs) for the latest state).

### Shortest Path Problems
* [Breadth-First Search](benchmarks/sssp/breadth_first_search)
* [Integral-Weight SSSP](benchmarks/sssp/integral_weight_sssp)
* [Positive-Weight SSSP](benchmarks/sssp/positive_weight_sssp)
* [General-Weight SSSP](benchmarks/sssp/general_weight_sssp)
* [Single-Source Widest Path](benchmarks/sssp/ss_widest_path)
* [Single-Source Betwenness Centrality](benchmarks/sssp/ss_betweenness_centrality)
* [Graph Spanner](benchmarks/sssp/spanner)

### Connectivity Problems
* [Low-Diameter Decomposition](benchmarks/connectivity/low_diameter_decomposition)
* [Connectivity](benchmarks/connectivity/connectivity)
* [Biconnectivity](benchmarks/connectivity/biconnectivity)
* [Strongly Connected Components](benchmarks/connectivity/strongly_connected_components)
* [Minimum Spanning Forest](benchmarks/connectivity/minimum_spanning_forest)

### Covering Problems
* [Maximal Independent Set (MIS)](benchmarks/covering/maximal_independent_set)
* [Maximal Matching](benchmarks/covering/maximal_matching)
* [Graph Coloring](benchmarks/covering/coloring)
* [Approximate Set Cover](benchmarks/covering/apx_set_cover)

### Substructure Problems
* [Triangle Counting](benchmarks/substructure/triangle_counting)
* [$k$-Clique Counting](benchmarks/substructure/k_clique_counting)
* [$k$-Core](benchmarks/substructure/k_core)
* [$k$-Truss](benchmarks/substructure/k_truss)
* [Approximate Densest Subgraph](benchmarks/substructure/apx_densest_subgraph)
* [Low-Outdegree Orientation (Arboricity Ordering)](benchmarks/substructure/low_outdegree_orientation)

### Eigenvector Problems
* [PageRank](benchmarks/eigenvector/pagerank)
* [CoSimRank](benchmarks/eigenvector/cosimrank)

