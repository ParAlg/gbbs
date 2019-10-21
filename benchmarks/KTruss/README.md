### Problem: KTruss

This directory contains a preliminary version of a k-Truss benchmark. The
algorithm is currently unpublished, but if you use this code in your work,
please consider citing our SPAA'18 paper.

The code uses bucketing, and works similarly to the k-core code available in
GBBS.

You can run the code as follows:

```
$ make
g++ -I../../ -mcx16 -ldl -std=c++17 -march=native -O3 -g -DLONG  -DAMORTIZEDPD  -DUSEMALLOC -DHOMEGROWN -pthread -o KTruss KTruss.cc

$ numactl -i all ./KTruss -m -nb 16 -rounds 1 -s ~/inputs/soc-LiveJournal1_sym.adj

$ numactl -i all ./KTruss -m -nb 16 -rounds 1 -s ~/inputs/twitter_sym.adj
```

* -s indicates that the input graph is symmetric
* -nb indicates the number of buckets to use (16 seems to be a reasonable choice
    for real-world graphs)
* -m mmaps the input file, which can help speed up graph loading

The trussness values are currently stored in an intermediate structure (array of
hashtables). We plan to update this benchmark in the coming month(s) with
optimizations and improvements, and will specify a better output format.

Note that due to an optimization that packs out neighbor-lists, the benchmark
only runs once. This restriction can be removed in general (e.g., by copying the
graph, or using more sophisticated techniques for filtering edges), but the
preliminary version currently does not support these features.

If you have questions, please contact {ldhulipa@cs.cmu.edu, jshun@mit.edu}.
