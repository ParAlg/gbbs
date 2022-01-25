# GBBS: Graph Based Benchmark Suite

Parallel Nucleus Decomposition
--------

This folder contains code for our parallel (r, s)-nucleus decomposition algorithms.
Detailed information about the required compilation system and 
input graph formats can be found in the top-level directory of this 
repository. We describe here the various options implemented from our 
paper, [Theoretically and Practically Efficient Parallel Nucleus Decomposition](https://arxiv.org/pdf/2111.10980.pdf).

Running Code
-------
First, note that for k-truss, or (2, 3)-nucleus decomposition, the most efficient
code is found in `gbbs/benchmarks/KTruss`. The graph contraction optimization
is also included in that directory.

The nucleus decomposition application take the input graph as input, as well as 
flags to specify desired optimizations. Note that the `-s` flag must be set to 
indicate a symmetric (undirected) graph.

The basic options for arguments are:
* `-r` followed by a positive integer > 2, which specifies r in (r, s)-nucleus
decomposition.
* `-ss` followed by a positive integer > r, which specifies s in (r, s)-nucleus
decomposition.
* `-relabel`, which if set indicates that the graph relabeling optimization
should be used.
* `-efficient` followed by `0`, `1`, or `2`, which specifies the method to
obtain the set of updated r-cliques. `0` corresponds to a simple array, `1`
corresponds to a list buffer, and `2` corresponds to a hash table.

Note that some optimizations are interdependent and only apply if other
optimizations are activated. Each of the following combinations of arguments can
be combined with any of the basic options (but not with each other), and refer
to the different possibilities for optimizations relating to the main parallel
hash table used in the algorithm. More details of these are discussed in the 
paper.

* `-tt 1`: One-level parallel hash table
* `-tt 2`: Two-level non-contiguous parallel hash table
* `-tt 3 -nl` followed by a positive integer `x` where `x`> 1: `x`-multi-level
non-contiguous parallel hash table 
* `-tt 2 -contig`: Two-level contiguous parallel hash table using binary search
* `-tt 3 -contig -nl` followed by a positive integer `x` where `x`> 1: `x`-multi-level
contiguous parallel hash table using binary search
* `-tt 5`: Two-level contiguous parallel hash table using stored pointers
* `-tt 4 -nl` followed by a positive integer `x` where `x`> 1: `x`-multi-level
contiguous parallel hash table using stored pointers

 **Example Usage**

The main executable is `NucleusDecomposition_main`. A template command is:

```sh
$ bazel run NucleusDecomposition_main -- -rounds 1 -s -r 3 -ss 4 -relabel -efficient 1 -tt 5 </path/to/input/graph>
```

<!---
Parallel Hash Table Options:

One-level

Two-level + Non-contiguous

\ell-mulit-level + Non-contiguous

Two-level + Contiguous + Binary search

\ell-multi-level + Contiguous + Binary search

Two-level + Contiguous + Stored pointers

\ell-multi-level + Contiguous + Stored pointers


Graph Relabeling:


Obtaining Set of Updated r-Cliques:

Simple Array

List Buffer

Hash Table


-tt 1 is one-level hash twotable
-tt 3 is multi-level hash table; set -nl to 2 or 3 (so set -nl to l where we have an l-multi-level)
-tt 2 is a two-level combination of an array and parallel hash table

the following options apply for -tt2 and -tt3 ( do not apply for -tt1 ):

for -tt2 and -tt3, you can set -contig for contiguous space, which automatically uses binary search

if you want to use contig with stored pointers, -tt2 should be -tt5, and -tt3 should be -tt4

-relabel for graph relabeling

-efficient can be set to 0 (simple array), 1 (list buffer), 2 (hash table)

-->