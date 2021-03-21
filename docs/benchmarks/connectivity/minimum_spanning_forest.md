---
id: minimum_spanning_forest
title: Minimum Spanning Forest
---


## Problem Specification
#### Input
$G=(V, E)$, an undirected, weighted graph on $n$ vertices.

#### Output
$T$, a set of edges representing a minimum spanning forest (MSF) of $G$.


## Algorithm Implementations
We present an implementation of Boruvka's algorithm that runs in
$O(m\log n)$ work and $O(\log^2 n)$ depth w.h.p.  A detailed
description of our algorithm can be found in Section 6.2 of [this
paper](https://ldhulipala.github.io/papers/gbbs_topc.pdf).

The code for the primary implementation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/MinimumSpanningForest/Boruvka).


## Cost Bounds

The algorithm runs in $O(m\log n)$ work and $O(\log^{2} n)$ depth
w.h.p.

