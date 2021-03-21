---
id: k_core
title: k-Core (Coreness)
---

## Problem Specification
#### Input
$G=(V, E)$, an undirected graph.

#### Output
$D$, a [mapping](/docs/benchmarks/definitions) from each vertex to its
coreness value.  See below for the definition of $k$-cores and coreness values.


## Definitions
A **$k$-core** of a graph is a maximal subgraph $H$ where the
degree of every vertex in $H$ is $\geq k$. The **coreness** of a
vertex is the maximum $k$-core a vertex participates in. The $k$-core
problem in this paper is to compute a mapping from each vertex to its
coreness value.


## Algorithm Implementations

We provide an implementation of $k$-core based on the work-efficient
peeling algorithm from Julienne [1] in GBBS. We provide a tutorial on
how to implement this $k$-core example in [our tutorial](tutorial/kcore_tutorial).

The code for our implemenation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/KCore/JulienneDBS17/).

## Cost Bounds

The algorithm runs in $O(n + m)$ work and $O(\rho(G) \log n)$
depth, where $\rho(G)$ is the peeling-complexity of the graph $G$,
defined as the number of rounds to peel the graph to an empty graph
where each peeling step removes all minimum degree vertices. More
details can be found in Section 6.4 of [2].


## Compiling and Running

The benchmark can be compiled by running:
```
bazel build -c opt //benchmarks/BFS/KCore/JulienneDBS17/...
```

It can be run on a test input graph in the *uncompressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/KCore/JulienneDBS17/KCore_main  -s -m -src 1 inputs/rMatGraph_J_5_100
```

It can be run on a test input graph in the *compressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/KCore/JulienneDBS17/KCore_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda
```

## References

[1] Laxman Dhulipala, Guy Blelloch, and Julian Shun<br/>
[*Julienne: A Framework for Parallel Graph Algorithms using Work-efficient Bucketing*](https://ldhulipala.github.io/papers/Bucketing.pdf)<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 293-304, 2017.

[2] Laxman Dhulipala, Guy Blelloch, and Julian Shun<br/>
[*Theoretically Efficient Parallel Graph Algorithms Can Be Fast and Scalable*](https://ldhulipala.github.io/papers/gbbs_topc.pdf)<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 393-404, 2018. <br/>
