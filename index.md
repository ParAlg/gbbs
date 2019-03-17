---
title: ""
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

To understand the performance and resource-requirements of graph
processing systems capable of processing very large graphs we maintain
a list of our times as well as times achieved by state-of-the-art
systems on a [Scoreboard](/scoreboard).  Please contact
<ldhulipa@cs.cmu.edu> to add your new experimental results to the
scoreboard. Our benchmark includes:

### Distance Problems
* Breadth-First Search
* Integral-Weight SSSP
* General-Weight SSSP
* Single-Source Betweenness Centrality

### Connectivity Problems
* Low-Diameter Decomposition
* Connectivity
* Biconnectivity
* Strongly Connected Components
* Minimum Spanning Forest

### Symmetry-Breaking Problems
* Maximal Matching
* Maximal Independent Set
* Graph Coloring

### Other Problems
* KCore
* Approximate Set Cover
* Triangle Counting

If you find our work useful, please consider citing us:

```
@inproceedings{dhulipala2018theoretically,
  author    = {Laxman Dhulipala and
               Guy E. Blelloch and
               Julian Shun},
  title     = {Theoretically Efficient Parallel Graph
               Algorithms Can Be Fast and Scalable},
  booktitle = {ACM Symposium on Parallelism in Algorithms
               and Architectures (SPAA)},
  year      = {2018},
}
```

