---
id: strongly_connected_components
title: Strongly Connected Components
---


## Problem Specification
#### Input
$G=(V, E)$, an undirected graph on $n$ vertices. The input graph can
either be weighted or unweighted.

#### Output
$\mathcal{L}$, a [mapping](/docs/benchmarks/definitions) from each
vertex to the label of its strongly connected component.

## Algorithm Implementations

We present the first implementation of the SCC algorithm from Blelloch
et al. [1].
The code for our implementation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/StronglyConnectedComponents/RandomGreedyBGSS16).


## Cost Bounds

The algorithm runs in $O(m\log n)$ expected work and
$O(\mathsf{diam}(G) \log n)$ depth w.h.p. A detailed proof of the
algorithm's theoretical costs can be found in [1].


## References

[1] Guy Blelloch, Yan Gu, Julian Shun, and Yihan Sun<br/>
*Parallelism in Randomized Incremental Algorithms*<br/>
Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 467-478, 2016.
