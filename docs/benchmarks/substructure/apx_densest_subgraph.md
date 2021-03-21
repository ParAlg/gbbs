---
id: apx_densest_subgraph
title: Approximate Densest Subgraph
---

## Problem Specification
#### Input
$G=(V, E)$, an undirected graph, and a parameter $\epsilon$.

#### Output
$U \subseteq V$, a set of vertices such that the density of $G_{U}$
is a $2(1+\epsilon)$ approximation of the density of the densest subgraph of $G$.


## Definitions

The **density** of a subset of vertices $S$ is the number of edges in
the subgraph $S$ divided by the number of vertices.

## Algorithm Implementations

In GBBS we implement the elegant $(2+\epsilon)$-approximation
algorithm of Bahmani et al. [1]. More details about our implementation
can be found in Section 6.4 of [2].

The code for our implemenation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/ApproximateDensestSubgraph/ApproxPeelingBKV12/).

## Cost Bounds

Our implementation of the algorithm runs in $O(m+n)$ work and
$O(\log_{1+\epsilon} n \cdot \log n)$ depth. We provide a proof of
this result in Section 6.4 of [2].


## Compiling and Running

The benchmark can be compiled by running:
```
bazel build -c opt //benchmarks/ApproximateDensestSubgraph/ApproxPeelingBKV12/...
```

It can be run on a test input graph in the *uncompressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/ApproximateDensestSubgraph/ApproxPeelingBKV12/DensestSubgraph_main -s -m -src 1 inputs/rMatGraph_J_5_100
```

It can be run on a test input graph in the *compressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/ApproximateDensestSubgraph/ApproxPeelingBKV12/DensestSubgraph_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda
```

## References

[1] Bahman Bahmani, Ravi Kumar, and Sergei Vassilvitskii<br/>
*Densest Subgraph in Streaming and MapReduce*<br/>
Proc. {VLDB} Endow. 5(5), pp. 454-465, 2012

[2] Laxman Dhulipala, Guy Blelloch, and Julian Shun<br/>
[*Theoretically Efficient Parallel Graph Algorithms Can Be Fast and Scalable*](https://ldhulipala.github.io/papers/gbbs_topc.pdf)<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 393-404, 2018. <br/>
