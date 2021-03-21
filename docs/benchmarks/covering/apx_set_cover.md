---
id: apx_set_cover
title: Approximate Set Cover
---

## Problem Specification
#### Input
$G=(V=(S,E), A)$, an undirected bipartite graph representing an unweighted set cover instance.

#### Output
$S' \subseteq S$, a set of sets such that $\cup_{s \in S'}
N(s) = E$, and $|S'|$ is an $O(\log n)$-approximation to the optimal cover.


## Algorithm Implementations

We provide an implementation of the Blelloch et al. algorithm [1]. Our
implementation is specifically based on the bucketing-based
implementation from Julienne [2], with one significant change
regarding how sets acquire elements which we discuss in our paper [3].

The code for our implemenation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/ApproximateSetCover/MANISBPT11).


## Cost Bounds
The algorithm runs in $O(m)$ work and $O(\log^3 n)$ depth. A detailed
proof can be found in [1], and a high-level description of the
algorithm can be found in Section 6.3 of [3].


## Compiling and Running

The benchmark can be compiled by running:
```
bazel build -c opt //benchmarks/ApproximateSetCover/MANISBPT11/...
```

It can then be run on a test input graph in the *uncompressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/ApproximateSetCover/MANISBPT11/ApproximateSetCover_main -s -m -src 1 inputs/rMatGraph_J_5_100
```

It can then be run on a test input graph in the *compressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/ApproximateSetCover/MANISBPT11/ApproximateSetCover_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda
```

## References

[1] Guy Blelloch, Richard Peng, and Kanat Tangwongsan<br/>
*Linear-work Greedy Parallel Approximate Set Cover and Variants*<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), 2011.

[2] Laxman Dhulipala, Guy Blelloch, and Julian Shun<br/>
*Julienne: A Framework for Parallel Graph Algorithms using Work-efficient Bucketing*<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 293-304, 2017.

[3] Laxman Dhulipala, Guy Blelloch, and Julian Shun<br/>
[*Theoretically Efficient Parallel Graph Algorithms Can Be Fast and Scalable*](https://ldhulipala.github.io/papers/gbbs_topc.pdf)<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 393-404, 2018. <br/>
