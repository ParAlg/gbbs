---
id: spanner
title: Graph Spanner
---

## Problem Specification
#### Input
$G=(V, E)$, an undirected, unweighted graph, and an integer
stretch factor, $k$.

#### Output
$H \subseteq E$, a set of edges such that for every $u,v \in V$
connected in $G$, $\mathsf{dist}_{H}(u, v) \leq O(k) \cdot \mathsf{dist}_{G}(u,v)$.

## Algorithm Implementations

The algorithm is based on the Miller, Peng, Xu, and Vladu (MPXV) paper from SPAA'15 [1].
Our implementation is described in more detail in [2].
The code for our implemenation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/Spanner/MPXV15).

The construction results in an $O(k)$-spanner with expected size
$O(n^{1+1/k})$.

## Cost Bounds

The algorithm runs in $O(m)$ work and $O(k\log n)$ depth w.h.p.
More details about our implementation and the cost bounds can be found in Section 6.1 of [2].


## Compiling and Running

The benchmark can be compiled by running:
```
bazel build -c opt //benchmarks/Spanner/MPXV15:Spanner
```

It can then be run on a test input graph in the *uncompressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/Spanner/MPXV15/Spanner_main -s -m -src 1 inputs/rMatGraph_J_5_100
```

It can then be run on a test input graph in the *compressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/Spanner/MPXV15/Spanner_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda
```

## References

[1] Gary L. Miller, Richard Peng, Adrian Vladu, Shen Chen Xu.<br/>
[*Improved Parallel Algorithms for Spanners and Hopsets*](https://arxiv.org/abs/1309.3545). <br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), 2015.

[2] Laxman Dhulipala, Guy Blelloch, and Julian Shun<br/>
[*Theoretically Efficient Parallel Graph Algorithms Can Be Fast and Scalable*](https://ldhulipala.github.io/papers/gbbs_topc.pdf)<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 393-404, 2018. <br/>
