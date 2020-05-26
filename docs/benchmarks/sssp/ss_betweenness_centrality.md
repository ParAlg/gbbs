---
id: ss_betweenness_centrality
title: Single-Source Betweenness Centrality
---

## Problem Specification
#### Input
$G=(V, E)$, an undirected graph, and a source $s \in V$. The input graph can either be undirected, or directed.

#### Output
Output: $S$, a mapping from each vertex $v$ to the centrality
contribution from all $(s, t)$ shortest paths that pass through $v$.
Please see below for the definition.

## Algorithm Implementations

Betweenness centrality is a classic tool in social network
analysis for measuring the importance of vertices in a
network.

We require several definitions to properly specify the benchmark. 
* Define
$\sigma_{st}(v)$ to be the total number of $s-t$ shortest paths,
$\sigma_{st}(v)$ to be the number of $s-t$ shortest paths that pass
through $v$, and $\delta_{st}(v) = \frac{\sigma_{st}(v)}{\sigma_{st}}$
to be the \defn{pair-dependency} of $s$ and $t$ on $v$. 
* The
*betweenness centrality* of a vertex $v$ is equal to $\sum_{s
\neq v \neq t \in V} \delta_{st}(v)$, i.e. the sum of
pair-dependencies of shortest-paths passing through $v$.
* Brandes proposes an algorithm to compute the
betweenness centrality values based on the dependency of a vertex: the
*dependency* of a vertex $r$ on a vertex $v$ is $\delta_{r}(v) = \sum_{t \in V} \delta_{rt}(v)$. 

The single-source betweenness
centrality benchmark is to compute these dependency values for each vertex given a source vertex.

The code for our implemenation is available
[here](https://github.com/ldhulipala/gbbs/tree/master/benchmarks/SSBetweennessCentrality/Brandes).

## Cost Bounds

The algorithm runs in $O(n + m)$ work and $O(\mathsf{Diam}(G) \log n)$
depth. Please [1] for details.


## Compiling and Running

The benchmark can be compiled by running:
```
bazel build -c opt //benchmarks/SSBetweennessCentrality/Brandes:SSBetweennessCentrality
```

It can then be run on a test input graph in the *uncompressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/SSBetweennessCentrality/Brandes/SSBetweennessCentrality_main -s -m -src 1 inputs/rMatGraph_J_5_100
```

It can then be run on a test input graph in the *compressed format* as follows:
```
numactl -i all ./bazel-bin/benchmarks/SSBetweennessCentrality/Brandes/SSBetweennessCentrality_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda
```

## References

[1] Laxman Dhulipala, Guy E. Blelloch, and Julian Shun. [Theoretically Efficient Parallel Graph Algorithms Can Be Fast and Scalable](https://arxiv.org/abs/1805.05208). Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 393-404, 2018.
