---
id: maximal_independent_set
title: Maximal Independent Set
---

## Problem Specification
#### Input
$G=(V, E)$, an undirected graph.

#### Output
$U \subseteq V$, a set of vertices such that no two vertices in $U$ are neighbors and all vertices in $V \setminus U$ have a neighbor in $U$.

## Algorithm Implementations
The maximal independent set problem is to compute a subset of vertices
$U$ such that no two vertices in $U$ are neighbors, and all vertices
in $V \setminus U$ have a neighbor in $U$.  Maximal independent set
(MIS) and maximal matching (MM) are easily solved in linear work
sequentially using greedy algorithms.

In GBBS, we implement the rootset-based algorithm for MIS from
Blelloch et al.[1] which runs in $O(m)$ work
and $O(\log^2 n)$ depth w.h.p. (using the
improved depth analysis of Fischer and Noever [2]). More details about
our algorithm implementation can be found in [3].

The code for our implemenation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/GraphColoring/Hasenplaugh14/).

## Cost Bounds
The algorithm runs in $O(m)$ work and $O(\log^2 n)$ depth w.h.p.

## Compiling and Running

The benchmark can be compiled by running:
```
bazel build -c opt //benchmarks/MaximalIndependentSet/RandomGreedy/...
```

It can then be run on a test input graph in the *uncompressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/MaximalIndependentSet/RandomGreedy/MaximalIndependentSet_main -s -m -src 1 inputs/rMatGraph_J_5_100
```

It can then be run on a test input graph in the *compressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/MaximalIndependentSet/RandomGreedy/MaximalIndependentSet_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda
```

## References

[1] Guy Blelloch, Jeremy Fineman, and Julian Shun<br/>
*Greedy Sequential Maximal Independent Set and Matching are Parallel on Average*<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 308-317, 2012.

[2] Manuella Fischer and Andreas Noever<br/>
*Tight Analysis of Parallel Randomized Greedy MIS*<br/>
Proceedings of the ACM-SIAM Symposium on Discrete Algorithms (SODA), pp.2152-2160, 2018

[3] Laxman Dhulipala, Guy Blelloch, and Julian Shun<br/>
[*Theoretically Efficient Parallel Graph Algorithms Can Be Fast and Scalable*](https://ldhulipala.github.io/papers/gbbs_topc.pdf)<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 393-404, 2018.
