---
id: spanning_forest
title: Spanning Forest
---


## Problem Specification
#### Input
$G=(V, E)$, an undirected graph on $n$ vertices.

#### Output
$T$, a set of edges representing a spanning forest of $G$.


## Algorithm Implementations
We provide multiple implementations of spanning forest in GBBS. The
primary implementation is based on the [low-diameter
decomposition](low_diameter_decomposition) based algorithm from Shun
et al. [1].

The code for the primary implementation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/SpanningForest/SDB14).


## Cost Bounds

The algorithm runs in $O(n + m)$ expected work and $O(\log^{3} n)$
depth w.h.p., and the proof can be found in the Shun et al. paper [1].


## Compiling and Running

The benchmark can be compiled by running:
```
bazel build -c opt //benchmarks/SpanningForest/SDB14/...
```

It can then be run on a test input graph in the *uncompressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/SpanningForest/SDB14/SpanningForest_main -s -m -src 1 inputs/rMatGraph_J_5_100
```

It can then be run on a test input graph in the *compressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/SpanningForest/SDB14/SpanningForest_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda
```

## References

[1] Julian Shun, Laxman Dhulipala, and Guy Blelloch<br/>
*A Simple and Practical Linear-Work Parallel Algorithm for Connectivity*<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 143-153, 2014.
