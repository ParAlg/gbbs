---
id: integral_weight_sssp
title: Integral-Weight SSSP (weighted BFS)
---

## Problem Specification
#### Input
$G=(V, E, w)$, a weighted graph with integral edge weights, and a
source, $s \in V$. The input graph can either be undirected, or
directed.

#### Output
Output: $D$, a mapping where $D[v]$ is the shortest path distance from
$s$ to $v$ in $G$ and $\infty$ if $v$ is unreachable.

## Algorithm Implementations
Our GBBS implementation implements the weighted breadth-first search
(wBFS) algorithm, a version of Dijkstra's algorithm that is well
suited for low-diameter graphs with small positive integer edge
weights. Our implementation uses the bucketing interface from
Julienne [1]. The idea of our algorithm is to maintain a bucket for
each possible distance, and to process them in order of increasing
distance. Each bucket is similar to a frontier in BFS, but unlike BFS,
when we process a neighbor $u$ of a vertex $v$ in the current bucket
$i$, we place $u$ in bucket $i + w_{uv}$.

The code for our implemenation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/IntegralWeightSSSP/JulienneDBS17).

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

[1] Laxman Dhulipala, Guy Blelloch, and Julian Shun<br/>
[*Julienne: A Framework for Parallel Graph Algorithms using Work-efficient Bucketing*](https://ldhulipala.github.io/papers/Bucketing.pdf)<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 293-304, 2017.

[2] Laxman Dhulipala, Guy Blelloch, and Julian Shun<br/>
[*Theoretically Efficient Parallel Graph Algorithms Can Be Fast and Scalable*](https://ldhulipala.github.io/papers/gbbs_topc.pdf)<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 393-404, 2018. <br/>
