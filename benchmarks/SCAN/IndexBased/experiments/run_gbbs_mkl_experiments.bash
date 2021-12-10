# Run GBBS parallel index SCAN experiments with MKL matrix multiplication

DATA_DIR=~/data
RESULTS_DIR=~/scan-experiment-results
RUN_PARALLEL="numactl -i all ./run_gbbs_mkl_experiments -s"
RUN_SERIAL="./run_gbbs_mkl_experiments -s"

mkdir --parents ${RESULTS_DIR}

declare -a UNWEIGHTED_GRAPHS=()
declare -a WEIGHTED_GRAPHS=("blood_vessel" "cochlea")

cd ..
make
cd -
for GRAPH in "${WEIGHTED_GRAPHS[@]}"
do
  ${RUN_PARALLEL} -wf ${DATA_DIR}/${GRAPH}-csr 2>&1 | tee ${RESULTS_DIR}/${GRAPH}-gbbs-mkl.txt
done

SERIAL=1
MKL_NUM_THREADS=1
cd ..
make
cd -
for GRAPH in "${WEIGHTED_GRAPHS[@]}"
do
  ${RUN_SERIAL} -wf ${DATA_DIR}/${GRAPH}-csr 2>&1 | tee ${RESULTS_DIR}/${GRAPH}-gbbs-mkl-serial.txt
done
