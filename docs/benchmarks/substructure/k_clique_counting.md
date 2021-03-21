---
id: k_clique_counting
title: kClique Counting
---

## Problem Specification
#### Input
$G=(V, E)$, an undirected graph.

#### Output
$K_{G}$, the total number of $k$-cliques in $G$. Each unordered $(v_1, \ldots, v_k)$ $k$-clique is counted once.

## Algorithm Implementations
In GBBS we implement the $k$-clique algorithm recently proposed by Shi
et al. [1]. The algorithm has polylogarithmic depth and is
work-efficient (matches the work of the best sequential algorithm) for
sparse graphs. The implementation is based on computing a
[low-outdegree orientation](low_outdegree_orientation).

The code for our implemenation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/CliqueCounting/).

## Cost Bounds

Our implementation runs in $O(m\cdot \alpha^{k-2})$ work and $O(\log^2
n)$ depth (for any constant $k$).
More details can be found in [1].


## Compiling and Running

The benchmark can be compiled by running:
```
bazel build -c opt //benchmarks/CliqueCounting/...
```

It can then be run on a test input graph in the *uncompressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/CliqueCounting/Clique_main -s -m -src 1 inputs/rMatGraph_J_5_100
```

It can then be run on a test input graph in the *compressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/CliqueCounting/Clique_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda
```

## References

[1] Jessica Shi, Laxman Dhulipala, and Julian Shun<br/>
*Parallel Clique Counting and Peeling Algorithms*<br/>
Under Submission<br/>
[arXiv Version](https://arxiv.org/abs/2002.10047)
