---
id: introduction
title: Introduction
description: Introduction to GBBS.
permalink: /docs/getting_started/introduction/
---

The Graph-Based Benchmark Suite (GBBS) is an active research project
to build *provably-efficient, scalable, single-machine* implementations
of *fundamental* graph algorithms. Our approach is to design *simple*
codes that rely on high-level functional APIs that have well-defined
cost-bounds in terms of the work-depth model. The project is in
development at both CMU and MIT, and has been used to implement over
20 problems so far.

## Getting Started
This part of the site provides instructions for installing, compiling, and
running the codes provided in GBBS. It also provide details on the graph
formats supported by GBBS.

* [Installation](install) provides instructions on
  system requirements and how to download and install GBBS.
* [Compilation](compile) explains how to compile benchmarks in
  GBBS using both our build-tool, `bazel`, and using `make`.
* [Graph Formats](formats) describes the graph formats
* supported by our implementations.
* [Running](run) explains how to run a benchmark on our
  graph inputs.
* [Graph Inputs](inputs) describes the main graph inputs used in our
   experiments, and provides details on how to obtain them.

## Benchmarks
This part of the site provides a problem-based benchmark
specifications for each of the problems supported by GBBS.
We provide an overview of the benchmarks in the [overview
page](/test_website/docs/benchmark/overview). Links to the benchmarks
can also be found below:

### Shortest Path Problems
* [Breadth-First Search](benchmarks/sssp/breadth_first_search)
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


## References

[1]
