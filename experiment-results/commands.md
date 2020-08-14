numactl -i all bazel run //benchmarks/SCAN/IndexBased:unweighted_experiments -- -s ~/data/orkut-gbbs 2>&1 | tee ~/gbbs/tomtseng/scan-experiments/experiment-results/orkut-gbbs.txt
bazel run --config=serial //benchmarks/SCAN/IndexBased:unweighted_experiments -- -s --serial ~/data/orkut-gbbs 2>&1 | tee ~/gbbs/tomtseng/scan-experiments/experiment-results/orkut-gbbs-serial.txt
