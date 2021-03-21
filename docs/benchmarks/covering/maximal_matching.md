---
id: maximal_matching
title: Maximal Matching
---

## Problem Specification
#### Input
$G=(V, E)$, an undirected graph.

#### Output
$E' \subseteq E$, a set of edges such that no two edges in $E'$ share an endpoint and all edges in $E \setminus E'$ share an endpoint with some edge in $E'$.

## Algorithm Implementations

In GBBS we implement a modified version of the prefix-based maximal
matching algorithm from Blelloch et al. [1].

We had to make several modifications to run the algorithm on the large
graphs in our experiments. The original code from [1] uses an edgelist
representation, but we cannot directly use this implementation as
uncompressing all edges would require a prohibitive amount of memory
for large graphs.  Instead, as in our MSF implementation, we simulate
the prefix-based approach by performing a constant number of
*filtering* steps. Each filter step packs out $3n/2$ of the highest
priority edges, randomly permutes them, and then runs the edgelist
based algorithm on the prefix.  After computing the new set of edges
that are added to the matching, we filter the remaining graph and
remove all edges that are incident to matched vertices. In practice,
just 3--4 filtering steps are sufficient to remove essentially all
edges in the graph. The last step uncompresses any remaining edges
into an edgelist and runs the prefix-based algorithm.  The filtering
steps can be done within the work and depth bounds of the original
algorithm.

The code for our implemenation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/BFS/NonDeterministicBFS).

## Cost Bounds

Our implementation of the Blelloch et al. algorithm runs in $O(m)$
expected work and $O(\log^2 m)$ depth w.h.p. (using the improved depth
shown in [2]).


## Compiling and Running

The benchmark can be compiled by running:
```
bazel build -c opt //benchmarks/MaximalMatching/RandomGreedy/...
```

It can then be run on a test input graph in the *uncompressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/MaximalMatching/RandomGreedy/MaximalMatching_main -s -m -src 1 inputs/rMatGraph_J_5_100
```

It can then be run on a test input graph in the *compressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/MaximalMatching/RandomGreedy/MaximalMatching_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda
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
