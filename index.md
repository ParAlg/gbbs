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
benchmark specifications and examples of how to run our codes and reproduce our
results.

We also maintain a [Scoreboard](/scoreboard/) which gathers running-times of
state-of-the-art graph processing systems capable of processing three of the
largest publicly available graphs. If you would like your experimental results
to be included please [contact us](mailto:ldhulipa@cs.cmu.edu,guyb@cs.cmu.edu,jshun@mit.edu).

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

If you find our work useful, please consider citing our paper:

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

