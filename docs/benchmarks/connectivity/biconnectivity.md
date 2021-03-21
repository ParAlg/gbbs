---
id: biconnectivity
title: Biconnectivity
---


## Problem Specification
#### Input
$G=(V, E)$, an undirected graph on $n$ vertices.

#### Output
$L$, a [mapping](/docs/benchmarks/definitions) from each edge to the
label of its biconnected component.


## Definitions

A *biconnected component* of an undirected graph is a maximal
subgraph such that the subgraph remains connected under the deletion
of any single vertex. Two closely related definitions are articulation
points and bridge. An *articulation point* is a vertex whose
deletion increases the number of connected components, and a
*bridge* is an edge whose deletion increases the number of
connected components.  Note that by definition an articulation point
must have degree greater than one.  The *biconnectivity problem* is to
emit a mapping that maps each edge to the label of its biconnected
component.

## Algorithm Implementations

We implement the Tarjan-Vishkin algorithm for biconnectivity. Our
implementation uses the parallel graph connectivity implementation
from GBBS, which is described
[here](/docs/benchmarks/connectivity/connectivity).

The code for the primary implementation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/Biconnectivity/TarjanVishkin).

## Cost Bounds

The algorithm runs in $O(n + m)$ expected work and
$O(\max(\mathsf{diam}(G)\log n, \log^3 n))$ depth w.h.p. A detailed
description of our algorithm can be found in Section 6.2 of [our
paper](https://ldhulipala.github.io/papers/gbbs_topc.pdf).

