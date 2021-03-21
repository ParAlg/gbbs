---
id: positive_weight_sssp
title: Positive-Weight SSSP (Delta Stepping)
---

## Problem Specification
#### Input
$G=(V, E, w)$, a weighted graph with positive edge weights, and a
source, $s \in V$. The input graph can either be undirected, or
directed.

#### Output
Output: $D$, a mapping where $D[v]$ is the shortest path distance from
$s$ to $v$ in $G$ and $\infty$ if $v$ is unreachable.

## Algorithm Implementations

The code for our implemenation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/PositiveWeightSSSP/DeltaStepping).
Our GBBS implementation is from the Julienne paper [1].

## Cost Bounds

Please [1] for details.


## Compiling and Running

The benchmark can be compiled by running:
```
bazel build -c opt //benchmarks/PositiveWeightSSSP/DeltaStepping:DeltaStepping
```

It can then be run on a test input graph in the *uncompressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/PositiveWeightSSSP/DeltaStepping/DeltaStepping_main -s -m -src 1 inputs/rMatGraph_J_5_100
```

It can then be run on a test input graph in the *compressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/PositiveWeightSSSP/DeltaStepping/DeltaStepping_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda
```

## References

[1] Laxman Dhulipala, Guy Blelloch, and Julian Shun<br/>
[*Julienne: A Framework for Parallel Graph Algorithms using Work-efficient Bucketing*](https://ldhulipala.github.io/papers/Bucketing.pdf)<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 293-304, 2017.

[2] Laxman Dhulipala, Guy Blelloch, and Julian Shun<br/>
[*Theoretically Efficient Parallel Graph Algorithms Can Be Fast and Scalable*](https://ldhulipala.github.io/papers/gbbs_topc.pdf)<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 393-404, 2018. <br/>
