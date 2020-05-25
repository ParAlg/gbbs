---
id: overview
title: Benchmarks Overview
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
repo](https://www.github.com/ldhulipala/gbbs) for the latest state).

### Shortest Path Problems
* [Breadth-First Search](sssp/breadth_first_search)
* [Integral-Weight SSSP](sssp/integral_weight_sssp)
* [Positive-Weight SSSP](sssp/positive_weight_sssp)
* [General-Weight SSSP](sssp/general_weight_sssp)
* [Single-Source Widest Path](sssp/ss_widest_path)
* [Single-Source Betwenness Centrality](sssp/ss_betweenness_centrality)
* [Graph Spanner](sssp/spanner)

### Connectivity Problems
* [Low-Diameter Decomposition](connectivity/low_diameter_decomposition)
* [Connectivity](connectivity/connectivity)
* [Biconnectivity](connectivity/biconnectivity)
* [Strongly Connected Components](connectivity/strongly_connected_components)
* [Minimum Spanning Forest](connectivity/minimum_spanning_forest)

### Covering Problems
* [Maximal Independent Set (MIS)](covering/maximal_independent_set)
* [Maximal Matching](covering/maximal_matching)
* [Graph Coloring](covering/coloring)
* [Approximate Set Cover](covering/apx_set_cover)

### Substructure Problems
* [Triangle Counting](substructure/triangle_counting)
* [$k$-Clique Counting](substructure/k_clique_counting)
* [$k$-Core](substructure/k_core)
* [$k$-Truss](substructure/k_truss)
* [Approximate Densest Subgraph](substructure/apx_densest_subgraph)
* [Low-Outdegree Orientation (Arboricity Ordering)](substructure/low_outdegree_orientation)

### Eigenvector Problems
* [PageRank](eigenvector/pagerank)
* [CoSimRank](eigenvector/cosimrank)


## References

[1]
