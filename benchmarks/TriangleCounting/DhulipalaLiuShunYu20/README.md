# Parallel Batch-Dynamic Triangle Counting
--------

This directory contains code for our experiments presented in the following paper: 

Laxman Dhulipala, Quanquan C. Liu, Julian Shun, Shangdi Yu: Parallel Batch-Dynamic k-Clique Counting. APOCS 2021: 129-143. [arXiv link](https://arxiv.org/abs/2003.13585)

Building
--------

Using make:
```
make
```

Using bazel:
```
bazel build //benchmarks/TriangleCounting/DhulipalaLiuShunYu20:DynamicTriangle_main
```
Input Graph Format
--------

GBBS uses the adjacency graph format based off [Problem Based Benchmark
suite](http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html)
and [Ligra](https://github.com/jshun/ligra).

For graphs obtained from [SNAP](https://snap.stanford.edu/snap/), you can use the [SNAP converter](https://github.com/ParAlg/gbbs/blob/tri_merge/utils/snap_converter.cc). Simply run the following commands from the [utils](https://github.com/ParAlg/gbbs/blob/tri_merge/utils/) directory:

```
# Get the graph from SNAP
wget https://snap.stanford.edu/data/wiki-Vote.txt.gz
gzip --decompress ${PWD}/wiki-Vote.txt.gz

# Build the converter
make snap_converter

# Use the converter to convert the desired file
./snap_converter -s -i <input file> -o <output file>
```
Note that the above converter can also be used to convert other graphs in edge list format.

Running the Code
--------

Command for running the experiments:
```
# dlsy insert
./Triangle -s -shuffle -w 1 -blocksize 128 -batchsize 1000000 -n 3072627  /ssd1/graphs/bench_experiments/com-orkut.ungraph.adj ../../../inputs/empty

# dlsy delete
./Triangle -s -shuffle -w 2 -blocksize 128 -batchsize 10000 -n 3072627  /ssd1/graphs/bench_experiments/com-orkut.ungraph.adj ../../../inputs/empty 

# makkar insert
./Triangle -s -shuffle -makkar -batchsize 100000 -n 3072627  /ssd1/graphs/bench_experiments/com-orkut.ungraph.adj ../../../inputs/empty 

 ```
 where "/ssd1/graphs/bench_experiments/com-orkut.ungraph.adj" should be replaced with the path to the input graph of choice in adjacency graph format (detailed above).
 
 Running RMAT:
 
 ```
 # DLSY DELETE
./Triangle -s -w 2 -blocksize 128 -batchsize 1000 -n 16384 ../../../inputs/empty  /ssd0/sy/rMatGraph_J_16K_1.6B
```

On NUMA machines, which can give improved performance:
```
# dlsy insert
numactl -i all ./Triangle -s -shuffle -w 1 -blocksize 128 -batchsize 1000000 -n 3072627  /ssd1/graphs/bench_experiments/com-orkut.ungraph.adj ../../../inputs/empty

# dlsy delete
numactl -i all ./Triangle -s -shuffle -w 2 -blocksize 128 -batchsize 10000 -n 3072627  /ssd1/graphs/bench_experiments/com-orkut.ungraph.adj ../../../inputs/empty 

# makkar insert
numactl -i all ./Triangle -s -shuffle -makkar -batchsize 100000 -n 3072627  /ssd1/graphs/bench_experiments/com-orkut.ungraph.adj ../../../inputs/empty 

 ```

 Running RMAT:
 ```
 # DLSY DELETE
numactl -i all ./Triangle -s -w 2 -blocksize 128 -batchsize 1000 -n 16384 ../../../inputs/empty  /ssd0/sy/rMatGraph_J_16K_1.6B
```

Running Twitter:
```
# dlsy delete
#numactl -i all ./Triangle -s -shuffle -w 2 -blocksize 128 -batchsize 1000  -n 41652231 /ssd1/graphs/bench_experiments/twitter_sym.adj ../../../inputs/empty 
```
