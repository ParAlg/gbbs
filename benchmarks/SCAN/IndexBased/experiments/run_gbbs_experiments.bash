# Run GBBS parallel index SCAN experiments.

DATA_DIR=~/data
RESULTS_DIR=~/scan-experiment-results
RUN_PARALLEL="numactl -i all bazel run --cxxopt=-O3 //benchmarks/SCAN/IndexBased/experiments:run_gbbs_experiments -- -s"
RUN_SERIAL="bazel run --cxxopt=-O3 --config=serial //benchmarks/SCAN/IndexBased/experiments:run_gbbs_experiments -- -s --serial"

mkdir --parents ${RESULTS_DIR}

declare -a UNWEIGHTED_GRAPHS=("orkut" "brain" "webbase" "friendster")
declare -a WEIGHTED_GRAPHS=("blood_vessel" "cochlea")

for GRAPH in "${UNWEIGHTED_GRAPHS[@]}"
do
  ${RUN_PARALLEL} ${DATA_DIR}/${GRAPH}-csr 2>&1 | tee ${RESULTS_DIR}/${GRAPH}-gbbs.txt
done
for GRAPH in "${WEIGHTED_GRAPHS[@]}"
do
  ${RUN_PARALLEL} -wf ${DATA_DIR}/${GRAPH}-csr 2>&1 | tee ${RESULTS_DIR}/${GRAPH}-gbbs.txt
done

for GRAPH in "${UNWEIGHTED_GRAPHS[@]}"
do
  ${RUN_SERIAL} ${DATA_DIR}/${GRAPH}-csr 2>&1 | tee ${RESULTS_DIR}/${GRAPH}-gbbs-serial.txt
done
for GRAPH in "${WEIGHTED_GRAPHS[@]}"
do
  ${RUN_SERIAL} -wf ${DATA_DIR}/${GRAPH}-csr 2>&1 | tee ${RESULTS_DIR}/${GRAPH}-gbbs-serial.txt
done
