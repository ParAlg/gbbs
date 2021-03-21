---
id: low_outdegree_orientation
title: Low-Outdegree Orientation
---

## Problem Specification
#### Input
$G=(V, E)$, an undirected graph.

#### Output
$S$, a total-ordering on the vertices s.t. orienting the edges of the
graph based on $S$ yields an $O(\alpha)-$orientation of the graph. See
below for the definitions of orientations.

## Definitions

The **arboricity**  ($\alpha$)} of a graph is the minimum number of
spanning forests needed to cover the graph. $\alpha$ is upper bounded
by $O(\sqrt{m})$ and lower bounded by $\Omega(1)$.

An **$c$-orientation** of an undirected graph is a total ordering on
the vertices, where the oriented out-degree of each vertex (the number
of its neighbors higher than it in the ordering) is bounded by $c$.

## Algorithm Implementations

We provide two implementations of low-outdegree orientation in GBBS. 
The implementations are based on the low-outdegree orientation
algorithms of Barenboim-Elkin and Goodrich-Pszona, which are efficient
in the $\mathsf{CONGEST}$ and I/O models of computation respectively.
We show that these algorithms can yield work-efficient and low-depth
algorithms in [3].

The code for our implemenations is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/DegeneracyOrder/).

## Cost Bounds

Both algorithms run in $O(n + m)$ work and $O(\log^2 n)$ depth w.h.p.
More details can be found in [3].


## References

[1] Leonid Barenboim and Michael Elkin<br/>
*Sublogarithmic Distributed MIS Algorithm for Sparse Graphs using Nash-Williams Decomposition*<br/>
Distributed Computing, 5(22), pp. 363-379, 2010

[2] Michael Goodrich and Pawel Pszona<br/>
*External-Memory Network Analysis Algorithms for Naturally Sparse Graphs*<br/>
Proceedings of the European Symposium on Algorithms, pp. 646-676, 2011

[3] Jessica Shi, Laxman Dhulipala, and Julian Shun<br/>
*Parallel Clique Counting and Peeling Algorithms*<br/>
Under Submission<br/>
[arXiv Version](https://arxiv.org/abs/2002.10047)
