---
id: general_weight_sssp
title: General-Weight SSSP (Bellman-Ford)
---

## Problem Specification
#### Input
$G=(V, E, w)$, a weighted graph, and a source, $s \in V$. The input graph can either be undirected, or directed.

#### Output
Output: $D$, a mapping where $D[v]$ is the shortest path distance from $s$ to $v$ in $G$ and $\infty$ if $v$ is unreachable.  If the graph contains any negative-weight cycles reachable from $s$, the vertices of these negative-weight cycles and vertices reachable from them must have a distance of $-\infty$.

## Algorithm Implementations

The code for our implemenation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/GeneralWeightSSSP/BellmanFord).
We provide more details about our implementation in [1].

## Cost Bounds

The algorithm runs in $O(\mathsf{Diam}(G)m)$ work and $O(\mathsf{Diam}(G) \log n)$
depth. Please [1] for details.


## Compiling and Running

The benchmark can be compiled by running:
```
bazel build -c opt //benchmarks/GeneralWeightSSSP/BellmanFord:BellmanFord_main
```

It can then be run on an input graph in the *uncompressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/GeneralWeightSSSP/BellmanFord/BellmanFord_main -s -m -src 1 inputs/rMatGraph_J_5_100
```

It can then be run on an input graph in the *compressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/GeneralWeightSSSP/BellmanFord/BellmanFord_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda
```

## References

[1] Laxman Dhulipala, Guy E. Blelloch, and Julian Shun. [Theoretically Efficient Parallel Graph Algorithms Can Be Fast and Scalable](https://arxiv.org/abs/1805.05208). Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 393-404, 2018.
