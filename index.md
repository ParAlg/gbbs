---
title: "The Graph Based Benchmark Suite (GBBS)"
layout: home
#header:
#  overlay_color: "#ccc"
#  overlay_filter: "0.5"
#  actions:
#    - label: "Download"
#      url: "https://github.com/ldhulipala/gbbs/"
---

This website supplements our paper,
["Theoretically Efficient Parallel Graph Algorithms Can Be Fast and
Scalable"](https://arxiv.org/abs/1805.05208), published in SPAA'18, with
benchmark specifications.

To encourage further research on efficient algorithms for very large graphs, we
have gathered performance numbers for graph processing systems capable of
solving these problems and present these times on our [Scoreboard](/scoreboard).
Please contact <ldhulipa@cs.cmu.edu> to add your experimental results to the
scoreboard.

Our benchmark consists of the following problems:

* Breadth-First Search
* Integral-Weight SSSP
* General-Weight SSSP
* Single-Source Betweenness Centrality
* Connectivity
* Biconnectivity
* Strongly Connected Components
* Minimum Spanning Forest
* Maximal Matching
* Maximal Independent Set
* Graph Coloring
* KCore
* Low-Diameter Decomposition
* Approximate Set Cover
* Triangle Counting

We also include specifications of the following extra graph problems:

* Approximate Densest Subgraph
* KTruss

If you find our work useful, please consider citing our
[paper](https://arxiv.org/abs/1805.05208):

```
@inproceedings{dhulipala2018theoretically,
  author    = {Laxman Dhulipala and
               Guy E. Blelloch and
               Julian Shun},
  title     = {Theoretically Efficient Parallel Graph Algorithms Can Be Fast and
               Scalable},
  booktitle = {ACM Symposium on Parallelism in Algorithms and Architectures (SPAA)},
  year      = {2018},
}
```

