---
id: definitions
title: Definitions
description: Definitions used in specifying GBBS Benchmarks
permalink: /docs/benchmarks/definitions
---

Here we specify the definitions of objects and terms used in our
benchmark specifications.

#### Unweighted Graphs
We consider both undirected (symmetric) and directed (asymmetric)
graphs in GBBS. We refer to graphs as $G=(V,E)$, where $V$ is a set of
vertices, and $E$ is a set of edges. Each edge consists of a pair of
vertices $(u,v)$.

When not specified, we assume that the underlying graph is *directed*.
The concrete formats supported by GBBS for unweighted graphs can be
found in the [formats page](/docs/formats)

#### Weighted Graphs
We also consider weighted graphs in GBBS, which we refer to as
$G=(V,E,w)$ where $w$ is a weighted function $w : E \mapsto
\mathcal{R}$. We refer to the weight of a $(u,v)$ edge as $w_{uv}$ The
benchmarks always specify the range, $\mathcal{R}$ of the weight
function, unless the algorithm is oblivious to it.

#### Distances
We use $\mathsf{dist}_{G}(u,v)$ to refer to the shortest path distance
between $u$ and $v$ in the graph $G$.

We use $\mathsf{diam}(G)$ to refer to the diameter of the graph, or the
longest shortest path distance in $G$ between any vertex $s$ and any
vertex $v$ reachable from $s$.


#### Mapping
Many of our benchmarks return $\alpha$-mapping, which are objects that
provide access to $\alpha$'s over a specified input domain. For
example, an `int`-mapping for the domain $[0, n)$ can be represented
* using an array of $n$ `int`s, or
* using a function from the domain to an `int` (implicitly represented).

The only requirement is that the returned object implements an array
operator, and can be accessed at every index in the domain.

