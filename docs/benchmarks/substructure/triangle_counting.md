---
id: triangle_counting
title: Triangle Counting
---

## Problem Specification
#### Input
$G=(V, E)$, an undirected graph.

#### Output
$T_{G}$, the total number of triangles in $G$. Each unordered $(u,v,w)$ triangle is counted once.

## Algorithm Implementations
In GBBS we implement the triangle counting algorithm described by Shun
and Tangwongsan [1].
This algorithm parallelizes Latapy's \emph{compact-forward}
algorithm [2], which creates a directed graph $DG$ where an edge
$(u,v) \in E$ is kept in $DG$ iff $\degree{u} < \degree{v}$.  Although
triangle counting can be done directly on the undirected graph in the
same work and depth asymptotically, directing the edges helps reduce
work, and ensures that every triangle is counted exactly once.
More details can be found in Section 6.4 of [3].

The code for our implemenation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/TriangleCounting/ShunTangwongsan15/).

## Cost Bounds

Our implementation runs in $O(m^{3/2})$ work and $O(\log n)$ depth.
More details can be found in Section 6.4 of [3].


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

[1] Julian Shun and Kanat Tangwongsan<br/>
*Multicore Triangle Computations Without Tuning*<br/>
Proceedings of the IEEE International Conference on Data Engineering (ICDE), pp. 149-160, 2015.

[2] Matthieu Latapy<br/>
*Main-memory Triangle Computations for Very Large (Sparse (Power-law)) Graphs*<br/>
Theor. Comput. Sci., 407(1-3), pp. 458-473, 2008

[3] Laxman Dhulipala, Guy Blelloch, and Julian Shun<br/>
[*Theoretically Efficient Parallel Graph Algorithms Can Be Fast and Scalable*](https://ldhulipala.github.io/papers/gbbs_topc.pdf)<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 393-404, 2018. <br/>
