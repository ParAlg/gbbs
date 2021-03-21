---
id: low_diameter_decomposition
title: Low Diameter Decomposition
---


## Problem Specification
#### Input
$G=(V, E)$, an undirected graph on $n$ vertices.

#### Output
$\mathcal{L}$, a [mapping](/docs/benchmarks/definitions) from each vertex to a cluster ID representing a $(O(\beta), O(\log n/\beta)$ decomposition.


## Definitions
A $(\beta, d)$-decomposition partitions $V$ into $C_{1}, \ldots, C_{k}$ such that:
* The shortest path between two vertices in $C_{i}$ using only vertices in $C_{i}$ is at most $d$.
* The number of edges $(u,v)$ where $u \in C_{i}, v \in C_{j}, j \neq i$ is at most $\beta m$.

The low-diameter decomposition problem that we study is to compute an
$(O(\beta), O(\log n/\beta)$ decomposition.


## Algorithm Implementations
We implement the Miller-Peng-Xu (MPX) algorithm [1] which computes an
$(2\beta, O(\log n / \beta))$ decomposition in $O(m)$ work and
$O(\log^2 n)$ depth w.h.p.  Our implementation is based on the
non-deterministic LDD implementation from Shun et al. [2] (designed as
part of a parallel connectivity implementation).

The code for the primary implementation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/LowDiameterDecomposition/MPX13).


## Cost Bounds

The algorithm runs in $O(n + m)$ expected work and $O(\log^{2} n)$
depth w.h.p. A detailed description of our algorithm can be found in
Section 6.2 of [this paper](https://ldhulipala.github.io/papers/gbbs_topc.pdf).


## References

[1] Gary L. Miller, Richard Peng, and Shen Chen Xu<br/>
*Parallel graph decompositions using random shifts*<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 196-203, 2013.

[2] Julian Shun, Laxman Dhulipala, and Guy Blelloch<br/>
*A Simple and Practical Linear-Work Parallel Algorithm for Connectivity*<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 143-153, 2014.

