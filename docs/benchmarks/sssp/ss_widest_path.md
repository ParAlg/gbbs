---
id: ss_widest_path
title: Single-Source Widest Path
---

## Problem Specification
#### Input
$G=(V, E, w)$, a weighted graph with integral edge weights, and a
source, $s \in V$. The input graph can either be undirected, or
directed.

#### Output
$D$, a mapping where $D[v]$ is the maximum over all paths
between $s$ and $v$ in $G$ of the minimum weight on the path and
$\infty$ if $v$ is unreachable.

## Algorithm Implementations

We provide two implementations of the algorithm, one based on the weighted SSWidestPath algorithm, and the other based on the Bellman-Ford algorithm.
We note that the Bellman-Ford implementation works for the version of this problem where the graph has general edge-weights.
The code for our implemenation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/SSWidestPath/JulienneDBS17).

## Cost Bounds

The wBFS-based implementation runs in $O(n + m)$ expected work and $O(\mathsf{Diam}(G) \log n)$ depth w.h.p. The Bellman-Ford based implementation runs in $O(\mathsf{Diam}(G)m)$ work and $O(\mathsf{Diam}(G) \log n)$ depth.


## Compiling and Running

The benchmark can be compiled by running:
```
bazel build -c opt //benchmarks/SSWidestPAth/JulienneDBS17:SSWidestPath
```

It can then be run on a test input graph in the *uncompressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/SSWidestPAth/JulienneDBS17/SSWidestPath_main -s -m -src 1 inputs/rMatGraph_J_5_100
```

It can then be run on a test input graph in the *compressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/SSWidestPAth/JulienneDBS17/SSWidestPath_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda
```

## References

[1] Laxman Dhulipala, Guy E. Blelloch, and Julian Shun. [Julienne: A Framework for Parallel Graph Algorithms using Work-efficient Bucketing](https://ldhulipala.github.io/papers/Bucketing.pdf). Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 293-304, 2017.

[2] Laxman Dhulipala, Guy E. Blelloch, and Julian Shun. [Theoretically Efficient Parallel Graph Algorithms Can Be Fast and Scalable](https://arxiv.org/abs/1805.05208). Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 393-404, 2018.
