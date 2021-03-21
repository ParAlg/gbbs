---
id: introduction
title: Introduction
description: Introduction to GBBS.
permalink: /docs/getting_started/introduction/
---

The Graph-Based Benchmark Suite (GBBS) is an active research project
to build *provably-efficient, scalable, single-machine* implementations
for a broad set of *fundamental* graph problems. Our approach is to
design *simple* algorithm implementations that rely on high-level
functional grpah and vertex primitives that have well-defined
cost-bounds in the work-depth model. The project is in development at
both CMU and MIT, and has been used to implement over 20 problems so
far.

## Getting Started
This part of the site provides instructions for installing, compiling, and
running the codes provided in GBBS. It also provide details on the graph
formats supported by GBBS.

* [Installation](install) provides instructions on
  system requirements and how to download and install GBBS.
* [Compilation](compile) explains how to compile benchmarks in
  GBBS using both our build-tool, `bazel`, and using `make`.
* [Graph Formats](formats) describes the graph formats
  supported by our implementations.
* [Running](run) explains how to run a benchmark on our
  graph inputs.
* [Graph Inputs](inputs) describes the main graph inputs used in our
   experiments, and provides details on how to obtain them.

## Benchmarks
This part of the site provides a problem-based benchmark
specifications for each of the problems supported by GBBS. Links to
the benchmarks can be found below:

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


## References

A list of publications built using GBBS and contributing to the
benchmark can be found [here](research).
