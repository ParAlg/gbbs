# GBBS: Graph Based Benchmark Suite

Clique Counting and Peeling Algorithms
--------

This folder contains code for our parallel k-clique counting and peeling.
Detailed information about the required compilation system and 
input graph formats can be found in the top-level directory of this 
repository. We describe here the various options implemented from our 
paper, [Parallel Clique Counting and Peeling Algorithms](https://arxiv.org/abs/2002.10047).


Running Code
-------
The applications take the input graph as input, as well as flags to specify
desired optimizations. Note that the `-s` flag must be set to indicate a symmetric 
(undirected) graph.

The options for arguments are:
* `-k` followed by a positive integer > 2, which specifies k.
* `--directType` followed by `GOODRICHPSZONA`, `BARENBOIMELKIN`, `KCORE`, `DEGREE`, or `ORIGINAL`, which 
specifies the ordering to use for the orientation in the algorithms. 
* `-e` followed by a float, which specifies the epsilon if the `GOODRICHPSZONA` or `BARENBOIMELKIN`
orderings are selected (otherwise, this argument is ignored).
* `--parallelType` followed by `VERT` or `EDGE`, which specifies if k-clique counting
should use parallelism per vertex or per edge.
* `--saveSpace`, which if set indicates that the first level of recursion should use
space proportional to the squared arboricity rather than linear space, to reduce space usage overall.
* `--peel`, which if set indicates that k-clique peeling should run following k-clique counting. Note 
that k-clique counting will compute counts per vertex, rather than in total.
* `--sparse`, which if set indicates that approximate k-clique counting should run. Note that
this cannot be set in conjunction with k-clique peeling or approximate k-clique peeling.
* `--colors` followed by a positive integer, which specifies the number of colors to 
use in approximate k-clique counting if `--sparse` is selected (otherwise, this argument is ignored).
* `--approxpeel`, which if set indicates that approximate k-clique peeling should run.
* `--approxeps` followed by a float, which specifies the epsilon to use for approximate k-clique
peeling if `--approxpeel` is selected (otherwise, this argument is ignored).

 **Example Usage**

The main executable is `Clique_main`. A template command is:

```sh
$ bazel run Clique_main -- -rounds 1 -k 4 --directType DEGREE --parallelType VERT --peel -s </path/to/input/graph>
```
