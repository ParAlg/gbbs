# GBBS: Graph Based Benchmark Suite

Five-Cycle Counting Algorithms
--------

This folder contains code for our parallel five-cycle counting algorithms.
Detailed information about the required compilation system and 
input graph formats can be found in the top-level directory of this 
repository. We describe here the various options implemented from our 
paper, [Parallel Five-Cycle Counting Algorithms](https://drops.dagstuhl.de/opus/volltexte/2021/13774/pdf/LIPIcs-SEA-2021-2.pdf).


Running Code
-------
The applications take the input graph as input, as well as flags to specify
desired optimizations. Note that the `-s` flag must be set to indicate a symmetric 
(undirected) graph.

The options for arguments are:
* `--escape`, which specifies that the ESCAPE-based algorithm should be run
  if set, and otherwise, the Kowalik-based algorithm is run.
* `--no-schedule`, which specifies that the work-scheduling optimization should
not be used if set, and otherwise, the work-scheduling optimization is used.
* `-o` followed by `0`, `1`, or `2`, which specifies the ordering to use for the arboricity orientation
in the algorithms. `0` corresponds to the parallel Goodrich-Pszona arboricity orientation, 
`1` corresponds to the parallel Barenboim-Elkin arboricity orientation, and `2`
corresponds to a degree orientation.
* `--serial`, which specifies that our serial implementation should be run
 if set, and otherwise, our parallel implementation is run.

 **Example Usage**

The main executable is `FiveCycle_main` in the `Parallel5Cycle/` directory.

After navigating to the `Parallel5Cycle/` directory, a template command is:

```sh
$ bazel run FiveCycle_main -- -rounds 1 -o 0 -s </path/to/input/graph>
```
