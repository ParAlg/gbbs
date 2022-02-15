# Parallel Batch-Dynamic Triangle Counting
--------

This directory contains code for our experiments presented in the following paper: 

Laxman Dhulipala, Quanquan C. Liu, Julian Shun, Shangdi Yu: Parallel Batch-Dynamic k-Clique Counting. APOCS 2021: 129-143. [arXiv link](https://arxiv.org/abs/2003.13585)

Building
--------

To build:
```
bazel build //benchmarks/TriangleCounting/DhulipalaLiuShunYu20:DynamicTriangle_main
```

Run Orkut:
```
# dlsy insert
numactl -i all ./Triangle -s -shuffle -w 1 -blocksize 128 -batchsize 1000000 -n 3072627  /ssd1/graphs/bench_experiments/com-orkut.ungraph.adj ../../../inputs/empty

# dlsy delelte
numactl -i all ./Triangle -s -shuffle -w 2 -blocksize 128 -batchsize 10000 -n 3072627  /ssd1/graphs/bench_experiments/com-orkut.ungraph.adj ../../../inputs/empty 

# makkar insert
numactl -i all ./Triangle -s -shuffle -makkar -batchsize 100000 -n 3072627  /ssd1/graphs/bench_experiments/com-orkut.ungraph.adj ../../../inputs/empty 

 ```

 Run RMAT:
 ```
 # DLSY DELETE
numactl -i all ./Triangle -s -w 2 -blocksize 128 -batchsize 1000 -n 16384 ../../../inputs/empty  /ssd0/sy/rMatGraph_J_16K_1.6B
```

Run Twitter:
```
# dlsy delete
#numactl -i all ./Triangle -s -shuffle -w 2 -blocksize 128 -batchsize 1000  -n 41652231 /ssd1/graphs/bench_experiments/twitter_sym.adj ../../../inputs/empty 
```
