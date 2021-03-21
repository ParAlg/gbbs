---
id: biconnectivity
title: Biconnectivity
---


## Problem Specification
#### Input
$G=(V, E)$, an undirected graph on $n$ vertices.

#### Output
$\mathcal{L}$, a [mapping](/docs/benchmarks/definitions) from each
edge to the label of its biconnected component.


## Definitions

A **biconnected component** of an undirected graph is a maximal
subgraph such that the subgraph remains connected under the deletion
of any single vertex. Two closely related definitions are articulation
points and bridge. An **articulation point** is a vertex whose
deletion increases the number of connected components, and a
**bridge** is an edge whose deletion increases the number of
connected components.  Note that by definition an articulation point
must have degree greater than one.  The **biconnectivity problem** is
to emit a mapping that maps each edge to the label of its biconnected
component.

## Algorithm Implementation

We implement the Tarjan-Vishkin algorithm for biconnectivity [1]. Our
implementation uses the parallel graph connectivity implementation as
a sub-routine. The connectivity benchmark is described
[here](/docs/benchmarks/connectivity/connectivity). We also make use
of the implicit representation of $\mathcal{L}$, the biconnectivity
labeling, proposed by Ben-David et al. [2].

The code for the primary implementation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/Biconnectivity/TarjanVishkin).

## Cost Bounds

The algorithm runs in $O(n + m)$ expected work and
$O(\max(\mathsf{diam}(G)\log n, \log^3 n))$ depth w.h.p. A detailed
description of our algorithm can be found in Section 6.2 of [this
paper](https://ldhulipala.github.io/papers/gbbs_topc.pdf).


## References

[1] Robert E. Tarjan and Uzi Vishkin<br/>
*An efficient parallel biconnectivity algorithm*<br/>
SIAM Journal on Computing, 14(4), pp. 862-874, 1985

[2] Naama Ben-David, Guy Blelloch, Jeremy Fineman, Phillip Gibbons, Yan Gu, Charles McGuffey, and Julian Shun<br/>
*Implicit Decomposition for Write-Efficient Connectivity Algorithms* <br/>
Proceedings of the IEEE International Parallel and Distributed Processing Symposium (IPDPS), pp. 711-722, 2018.
