---
id: pagerank
title: PageRank
---

## Problem Specification
#### Input
$G=(V, E)$, an undirected graph.

#### Output
$\mathcal{P}$, a [mapping](/docs/benchmarks/definitions) from each
vertex to its PageRank value after a single iteration of PageRank.

## Algorithm Implementations

Our PageRank algorithm terminates once the $l_{1}$ norm of the
differences between PageRank values between iterations is below
$\epsilon$. The algorithm that we implement is an extension of
the implementation of PageRank described in Ligra. More details about
our implementation can be found in Section 6.5 of [1].

The code for our implemenation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/PageRank/).

## Cost Bounds

The algorithm runs in $O(n + m)$ work and $O(\log n)$ depth for a
single iteration of PageRank.


## Compiling and Running

The benchmark can be compiled by running:
```
bazel build -c opt //benchmarks/PageRank/...
```

It can then be run on a test input graph in the *uncompressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/PageRank/PageRank_main -s -m -src 1 inputs/rMatGraph_J_5_100
```

It can then be run on a test input graph in the *compressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/PageRank/PageRank_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda
```

## References

[1] Laxman Dhulipala, Guy Blelloch, and Julian Shun<br/>
[*Theoretically Efficient Parallel Graph Algorithms Can Be Fast and Scalable*](https://ldhulipala.github.io/papers/gbbs_topc.pdf)<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 393-404, 2018. <br/>
