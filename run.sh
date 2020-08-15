RESULTS=~/gbbs/tomtseng/scan-experiments/experiment-results

#######################
# ppSCAN
#######################
cd ~/ppSCAN/pSCAN-refactor/build

# mkdir ~/data/orkut-pp && converter/converter ~/data/orkut.edges ~/data/orkut-pp/b_degree.bin ~/data/orkut-pp/b_adj.bin
numactl -i all ./pSCANParallelAVX2 ~/data/orkut-pp/ 0.6 5 2>&1 | tee ${RESULTS}/orkut-pp.txt
numactl -i all ./pSCANParallelAVX2 ~/data/brain-pp/ 0.6 5 2>&1 | tee ${RESULTS}/brain-pp.txt
mkdir ~/data/webbase-pp && converter/converter ~/data/webbase.edges ~/data/webbase-pp/b_degree.bin ~/data/webbase-pp/b_adj.bin
numactl -i all ./pSCANParallelAVX2 ~/data/webbase-pp/ 0.6 5 2>&1 | tee ${RESULTS}/webbase-pp.txt
mkdir ~/data/friendster-pp && converter/converter ~/data/friendster.edges ~/data/friendster-pp/b_degree.bin ~/data/friendster-pp/b_adj.bin
numactl -i all ./pSCANParallelAVX2 ~/data/friendster-pp/ 0.6 5 2>&1 | tee ${RESULTS}/friendster-pp.txt
# no blood_vessel or cochlea because those are weighted

##############################
# sequential index-based pt 1
###############################
cd ~/index-based-scan

mkdir -p ~/data/orkut-index
./run -f ~/data/orkut-index/ ~/data/orkut.edges
for i in {1..4}
do
  ./run -idx1 ~/data/orkut-index/ 2>&1 | tee --append ${RESULTS}/orkut-index.txt
  rm ~/data/orkut-index/idx_1
done
./run -idx1 ~/data/orkut-index/ -q 0.6 5 2>&1 | tee --append ${RESULTS}/orkut-index.txt
# skipping the other graphs for now since they'll take a long time

#######################
# GBBS experiments pt 1
#######################
cd ~/gbbs/tomtseng/scan-experiments

# redo timings with -O3. the --serial flag is on purpose --- it skips all the non-timing stuff
numactl -i all bazel run //benchmarks/SCAN/IndexBased:unweighted_experiments -- -s --serial ~/data/orkut-gbbs 2>&1 | tee ${RESULTS}/orkut-gbbs-2.txt
numactl -i all bazel run //benchmarks/SCAN/IndexBased:unweighted_experiments -- -s --serial ~/data/brain-gbbs 2>&1 | tee ${RESULTS}/brain-gbbs-2.txt
numactl -i all bazel run //benchmarks/SCAN/IndexBased:weighted_experiments -- -s -w --serial ~/data/blood_vessel-gbbs 2>&1 | tee ${RESULTS}/blood_vessel-gbbs-2.txt
numactl -i all bazel run //benchmarks/SCAN/IndexBased:weighted_experiments -- -s -w --serial ~/data/cochlea-gbbs 2>&1 | tee ${RESULTS}/cochlea-gbbs-2.txt

# webbase, friendster timings from scratch
numactl -i all bazel run //benchmarks/SCAN/IndexBased:unweighted_experiments -- -s ~/data/webbase-gbbs 2>&1 | tee ${RESULTS}/webbase-gbbs.txt
numactl -i all bazel run //benchmarks/SCAN/IndexBased:unweighted_experiments -- -s ~/data/friendster-gbbs 2>&1 | tee ${RESULTS}/friendster-gbbs.txt

# redo timing on orkut serial too
numactl -i all bazel run --config=serial //benchmarks/SCAN/IndexBased:unweighted_experiments -- -s --serial ~/data/orkut-gbbs 2>&1 | tee ${RESULTS}/orkut-gbbs-serial-2.txt
# let's not do other serial stuff yet --- it'll take too long

##############################
# sequential index-based pt 2
###############################
cd ~/index-based-scan

mkdir -p ~/data/brain-index
./run -f ~/data/brain-index/ ~/data/brain.edges
for i in {1..4}
do
  ./run -idx1 ~/data/brain-index/ 2>&1 | tee --append ${RESULTS}/brain-index.txt
  rm ~/data/brain-index/idx_1
done
./run -idx1 ~/data/brain-index/ -q 0.6 5 2>&1 | tee --append ${RESULTS}/brain-index.txt

mkdir -p ~/data/webbase-index
./run -f ~/data/webbase-index/ ~/data/webbase.edges
for i in {1..4}
do
  ./run -idx1 ~/data/webbase-index/ 2>&1 | tee --append ${RESULTS}/webbase-index.txt
  rm ~/data/webbase-index/idx_1
done
./run -idx1 ~/data/webbase-index/ -q 0.6 5 2>&1 | tee --append ${RESULTS}/webbase-index.txt

mkdir -p ~/data/friendster-index
./run -f ~/data/friendster-index/ ~/data/friendster.edges
for i in {1..4}
do
  ./run -idx1 ~/data/friendster-index/ 2>&1 | tee --append ${RESULTS}/friendster-index.txt
  rm ~/data/friendster-index/idx_1
done
./run -idx1 ~/data/friendster-index/ -q 0.6 5 2>&1 | tee --append ${RESULTS}/friendster-index.txt

#######################
# GBBS experiments pt 2
#######################
cd ~/gbbs/tomtseng/scan-experiments

# other serial stuff
numactl -i all bazel run --config=serial //benchmarks/SCAN/IndexBased:unweighted_experiments -- -s --serial ~/data/webbase-gbbs 2>&1 | tee ${RESULTS}/webbase-gbbs-serial.txt
numactl -i all bazel run --config=serial //benchmarks/SCAN/IndexBased:unweighted_experiments -- -s --serial ~/data/brain-gbbs 2>&1 | tee ${RESULTS}/brain-gbbs-serial.txt
numactl -i all bazel run --config=serial //benchmarks/SCAN/IndexBased:unweighted_experiments -- -s --serial ~/data/friendster-gbbs 2>&1 | tee ${RESULTS}/friendster-gbbs-serial.txt
numactl -i all bazel run --config=serial //benchmarks/SCAN/IndexBased:weighted_experiments -- -s -w --serial ~/data/blood_vessel-gbbs 2>&1 | tee ${RESULTS}/blood_vessel-gbbs-serial.txt
numactl -i all bazel run --config=serial //benchmarks/SCAN/IndexBased:weighted_experiments -- -s -w --serial ~/data/cochlea-gbbs 2>&1 | tee ${RESULTS}/cochlea-gbbs-serial.txt
