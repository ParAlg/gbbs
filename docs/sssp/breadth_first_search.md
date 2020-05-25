---
id: breadth_first_search
title: BFS
---

## Problem Specification
#### Input
$G=(V, E)$, an unweighted graph, and a source, $s \in V$. The input
graph can either be undirected, or directed.

#### Output
$P$, a [mapping](/benchmarks/definitions/) where $P[v]$ is the parent
of $v$ in the output BFS-tree rooted at $s$, and $P[s] = s$.

## Algorithm Implementations

We provide a single implementation of BFS in GAB. The implementation
is based on the non-deterministic BFS implementation in GAB. We
provide a tutorial on how to implement this BFS example in [our
tutorial](tutorial/bfs_tutorial).

The code for our implemenation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/BFS/NonDeterministicBFS).

## Cost Bounds

The algorithm runs in $O(n + m)$ work and $O(\mathsf{D}(G) \log n)$
depth, which follows from our bounds on the `edgeMap` primitive.


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

