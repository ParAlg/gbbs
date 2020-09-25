This directory contains code used for running experiments that were part of a
paper submission. The experiments measure the running time and the clustering
quality of the parallel index-based SCAN algorithm, the running time of
[ppSCAN](https://github.com/RapidsAtHKUST/ppSCAN),
and the running time of [GS\*-Index](http://www.vldb.org/pvldb/vol11/p243-wen.pdf).

Times from the original experiment runs used in the paper submission are in
`experiment-results.tar.gz`, which may be opened with the command `tar xvf
experiment-results.tar.gz`

To run the experiments:
* Run `bash download_graphs.bash` to download all the graphs used in the
  experiments.
* Run `bash run_gbbs_experiments.bash` to run the experiments with the GBBS
  index-based SCAN implementation.
* Run `bash run_ppscan_experiments.bash` to run timing experiments with ppSCAN.

We omit the code for running timing experiments with GS\*-Index since the
GS\*-Index code is not publicly available. We received the source code for
GS\*-Index via personal correspondence with its authors.
