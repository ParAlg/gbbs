CODE_DIR=~/scan-experiment-results/code
COMMIT_SHA=691b39309da1c6b5df46b264b5a300a35d644f70
DATA_DIR=~/data
GBBS_DIR=$(git rev-parse --show-toplevel)
RESULTS_DIR=~/scan-experiment-results

mkdir --parents ${CODE_DIR}
mkdir --parents ${RESULTS_DIR}

cd ${CODE_DIR}
git clone https://github.com/RapidsAtHKUST/ppSCAN.git
cd ppSCAN
git checkout ${COMMIT_SHA}
cd pSCAN-refactor
tar --extract --verbose --file=${GBBS_DIR}/benchmarks/SCAN/IndexBased/experiments/ppscan-changed-files.tar.gz
cp ppscan-changed-files/*.{h,cpp} .
mkdir -p build
cd build
cmake .. -DAVX2=ON -DKNL=OFF
make

declare -a GRAPHS=("orkut" "brain" "webbase" "friendster")

for GRAPH in "${GRAPHS[@]}"
do
  mkdir --parents ${DATA_DIR}/${GRAPH}-pp
  converter/converter ${DATA_DIR}/${GRAPH}-edges ${DATA_DIR}/${GRAPH}-pp/b_degree.bin ${DATA_DIR}/${GRAPH}-pp/b_adj.bin
  numactl -i all ./pSCANParallelAVX2 ${DATA_DIR}/${GRAPH}-pp/ 0.6 5 2>&1 >/dev/null | tee ${RESULTS_DIR}/${GRAPH}-pp.txt
done
