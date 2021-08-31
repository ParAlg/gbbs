# ConnectIt

This readme will contain instructions for how to use the ConnectIt framework
(VLDB'21).

1. Build the binaries using bazel: `bazel build //benchmarks/Connectivity/...`
2. Run a specific benchmark using one of the generated binaries. Some examples illustrating the usage follow. Note that the bazel-bin directory is located in the **top-level GBBS directory**.
  *  Run the unite_rem_cas algorithm without sampling: `numactl -i all ./bazel-bin/benchmarks/Connectivity/ConnectIt/mains/unite_rem_cas_nosample -rounds 3 -s  -m -r 10 -src 10012 ~/inputs/twitter_sym.adj`
  *  Run the unite_rem_cas algorithm with k-out sampling: `numactl -i all ./bazel-bin/benchmarks/Connectivity/ConnectIt/mains/unite_rem_cas_kout -rounds 3 -s  -m -r 10 -src 10012 ~/inputs/twitter_sym.adj`
  *  Run the unite_rem_cas algorithm with BFS sampling: `numactl -i all ./bazel-bin/benchmarks/Connectivity/ConnectIt/mains/unite_rem_cas_bfs -rounds 3 -s  -m -r 10 -src 10012 ~/inputs/twitter_sym.adj`
  *  Run the unite_rem_cas algorithm with LDD sampling: `numactl -i all ./bazel-bin/benchmarks/Connectivity/ConnectIt/mains/unite_rem_cas_ldd -rounds 3 -s  -m -r 10 -src 10012 ~/inputs/twitter_sym.adj`

To run on compressed graph inputs you should specify the `-c` flag. The -m flag uses mmap to save some time when reading the graph. Please see the main GBBS readme for more info.
