---
id: graph_coloring
title: Graph Coloring
---

## Problem Specification
#### Input
$G=(V, E)$, an undirected graph.

#### Output
$C$, a [mapping](/docs/benchmarks/definitions) from each vertex to a color such that for each
edge $(u, v) \in E$, $C(u) \neq C(v)$, using at most $\Delta+1$ colors.

## Algorithm Implementations

The graph coloring problem is to compute a mapping from each $v \in V$
to a color such that for each edge $(u, v) \in E$, $C(u) \neq C(v)$,
using at most $\Delta+1$ colors.  As graph coloring is
$\mathsf{NP}$-hard to solve optimally, algorithms like greedy
coloring, which guarantees a $(\Delta+1)$-coloring, are used instead
in practice, and often use much fewer than $(\Delta + 1)$ colors on
real-world graphs.

Jones and Plassmann (JP) [1] parallelize the greedy algorithm using linear
work, but unfortunately adversarial inputs exist for the heuristics
they consider that may force the algorithm to run in $O(n)$ depth.
Hasenplaugh et al. [2] introduce several heuristics that produce
high-quality colorings in practice and also achieve provably low-depth
regardless of the input graph. These include:
* LLF (largest-log-degree-first), which processes vertices ordered by the
log of their degree and
* SLL (smallest-log-degree-last), which processes vertices by removing all lowest log-degree vertices from the
graph, coloring the remaining graph, and finally coloring the removed
vertices.

We implement a synchronous version of the Jones-Plassmann algorithm
using the LLF heuristic. More details about our algorithm can be found
in Section 6.3 of [3].

The code for our implemenation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/GraphColoring/Hasenplaugh14/).

## Cost Bounds
The algorithm runs in $O(m + n)$ work and $O(L \log
\Delta + \log n)$ dept,h  where $L=\min\{\sqrt{m}, \Delta\} +
\log^{2}\Delta \log{n}/\log\log{n}$ in expectation [2,3].

## Compiling and Running

The benchmark can be compiled by running:
```
bazel build -c opt //benchmarks/GraphColoring/Hasenplaugh14/...
```

It can then be run on a test input graph in the *uncompressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/GraphColoring/Hasenplaugh14/GraphColoring_main -s -m -src 1 inputs/rMatGraph_J_5_100
```

It can then be run on a test input graph in the *compressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/GraphColoring/Hasenplaugh14/GraphColoring_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda
```

## References
[1] Mark Jones and Paul Plassmann<br/>
*A Parallel Graph Coloring Heuristic*<br/>
SIAM Journal on Scientific Computing, 1993

[2] William Hasenplaugh, Tim Kaler, Tao B. Schardl, and Charles Leiserson<br/>
*Ordering Heuristics for Parallel Graph Coloring*<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 166-177, 2014. <br/>

[3] Laxman Dhulipala, Guy Blelloch, and Julian Shun<br/>
[*Theoretically Efficient Parallel Graph Algorithms Can Be Fast and Scalable*](https://ldhulipala.github.io/papers/gbbs_topc.pdf)<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 393-404, 2018. <br/>
