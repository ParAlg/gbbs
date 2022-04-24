This directory contains code used for running experiments for the SIGMOD 2021
paper _Parallel Index-Based Structural Graph Clustering and Its Approximation_.
The experiments measure the running time and the clustering quality of the
parallel index-based SCAN algorithm, the running time of
[ppSCAN](https://github.com/RapidsAtHKUST/ppSCAN), and the running time of
[GS\*-Index](http://www.vldb.org/pvldb/vol11/p243-wen.pdf).

Times from the original experiment runs used in the paper submission are in
`experiment-results.tar.gz`, which may be opened with the command `tar xvf
experiment-results.tar.gz`

Beware: the experiments measuring the clustering quality are unoptimized and
slow. The experiments measuring serial running time are also extremely slow,
e.g., taking dozens of hours for the `cochlea` graph.

Installing dependencies:
* `sudo apt install python3.6-dev python3-pip --yes`
* `pip3 install sklearn`

To run the experiments:
* Run `bash download_graphs.bash` to download all the graphs used in the
  experiments.
* Run `bash run_gbbs_experiments.bash` to run the experiments with the GBBS
  index-based SCAN implementation.
* Run `bash run_ppscan_experiments.bash` to run timing experiments with ppSCAN.
* To run experiments using matrix multiplication, see the `tomtseng/mkl-scan`
  branch.

We omit the code for running timing experiments with GS\*-Index since the
GS\*-Index code is not publicly available. We received the source code for
GS\*-Index via personal correspondence with its authors.

For the paper, we ran the timing experiments on a c5.24xlarge AWS EC2 instance,
but we ran the clustering quality experiments on a m5a.24xlarge AWS EC2 instance
because the program barely exceeds the memory limits of a c5.24xlarge instance
when running on the [Friendster](https://snap.stanford.edu/data/com-Friendster.html) graph.
