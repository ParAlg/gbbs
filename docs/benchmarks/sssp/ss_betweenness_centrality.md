---
id: ss_betweenness_centrality
title: Single-Source Betweenness Centrality
---

## Problem Specification
#### Input
$G=(V, E, w)$, a weighted graph with positive edge weights, and a
source, $s \in V$. The input graph can either be undirected, or
directed.

#### Output
Output: $D$, a mapping where $D[v]$ is the shortest path distance from
$s$ to $v$ in $G$ and $\infty$ if $v$ is unreachable.

## Algorithm Implementations

The code for our implemenation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/PositiveWeightSSSP/DeltaStepping).

## Cost Bounds

The algorithm runs in $O(n + m)$ work and $O(\mathsf{Diam}(G) \log n)$
depth. Please [1] and [2] for details.


## Compiling and Running

The benchmark can be compiled by running:
```
bazel build -c opt //benchmarks/BFS/NonDeterministicBFS:BFS
```

It can then be run on a test input graph in the *uncompressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/BFS/NonDeterministicBFS/BFS_main -s -m -src 1 inputs/rMatGraph_J_5_100
```

It can then be run on a test input graph in the *compressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/BFS/NonDeterministicBFS/BFS_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda
```

## References

[1] Laxman Dhulipala, Guy E. Blelloch, and Julian Shun. [Julienne: A Framework for Parallel Graph Algorithms using Work-efficient Bucketing](https://ldhulipala.github.io/papers/Bucketing.pdf). Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 293-304, 2017.

[2] Laxman Dhulipala, Guy E. Blelloch, and Julian Shun. [Theoretically Efficient Parallel Graph Algorithms Can Be Fast and Scalable](https://arxiv.org/abs/1805.05208). Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 393-404, 2018.
