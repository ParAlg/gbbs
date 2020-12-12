# Download and format experimental graphs.

DATA_DIR=~/data
GBBS_DIR=$(git rev-parse --show-toplevel)

mkdir --parents ${DATA_DIR}
cd ${DATA_DIR}

echo "[Download Orkut graph]"
wget https://snap.stanford.edu/data/bigdata/communities/com-orkut.ungraph.txt.gz
gunzip com-orkut.ungraph.txt.gz
mv com-orkut.ungraph.txt orkut-edges

echo "[Download Friendster graph]"
wget https://snap.stanford.edu/data/bigdata/communities/com-friendster.ungraph.txt.gz
gunzip com-friendster.ungraph.txt.gz
mv com-friendster.ungraph.txt friendster-edges-uncompact

echo "[Download brain graph]"
wget http://nrvis.com/download/data/bn/bn-human-Jung2015_M87113878.zip
unzip bn-human-Jung2015_M87113878.zip
rm readme.html bn-human-Jung2015_M87113878.zip
mv bn-human-Jung2015_M87113878.edges brain-edges-uncompact

echo "[Download WebBase graph]"
wget https://www.cise.ufl.edu/research/sparse/MM/LAW/webbase-2001.tar.gz
tar --extract --verbose --file=webbase-2001.tar.gz
rm webbase-2001.tar.gz
mv webbase-2001/webbase-2001.mtx webbase-edges
# Remove the first 54 lines of the webbase edge list. The webbase file has 53
# lines of comments and another line giving the number of vertices and edges.
# These lines are unnecessary.
sed -i -e 1,54d webbase-edges
rmdir webbase-2001/

echo "[Download blood vessel graph]"
wget https://s3-us-west-2.amazonaws.com/humanbase/networks/blood_vessel_top.gz
gunzip blood_vessel_top.gz
mv blood_vessel_top blood_vessel-edges-uncompact

echo "[Download cochlea graph]"
wget https://s3-us-west-2.amazonaws.com/humanbase/networks/cochlea_top.gz
gunzip cochlea_top.gz
mv cochlea_top cochlea-edges-uncompact

cd ${GBBS_DIR}

EDGE_TO_CSR="numactl -i all bazel run //utils:snap_converter -- -s"
REMOVE_SINGLETONS="bazel run //benchmarks/SCAN/IndexBased/experiments:remove_singleton_vertices --"

echo "[Convert Orkut graph]"
${EDGE_TO_CSR} -i ${DATA_DIR}/orkut-edges -o ${DATA_DIR}/orkut-csr

echo "[Convert Friendster graph]"
${REMOVE_SINGLETONS} --skip 4 -i ${DATA_DIR}/friendster-edges-uncompact -o ${DATA_DIR}/friendster-edges
rm ${DATA_DIR}/friendster-edges-uncompact
${EDGE_TO_CSR} -i ${DATA_DIR}/friendster-edges -o ${DATA_DIR}/friendster-csr

echo "[Convert brain graph]"
${REMOVE_SINGLETONS} -i ${DATA_DIR}/brain-edges-uncompact -o ${DATA_DIR}/brain-edges
rm ${DATA_DIR}/brain-edges-uncompact
${EDGE_TO_CSR} -i ${DATA_DIR}/brain-edges -o ${DATA_DIR}/brain-csr

echo "[Convert WebBase graph]"
${EDGE_TO_CSR} -i ${DATA_DIR}/webbase-edges -o ${DATA_DIR}/webbase-csr

echo "[Convert blood vessel graph]"
${REMOVE_SINGLETONS} -wf -i ${DATA_DIR}/blood_vessel-edges-uncompact -o ${DATA_DIR}/blood_vessel-edges
rm ${DATA_DIR}/blood_vessel-edges-uncompact
${EDGE_TO_CSR} -wf -i ${DATA_DIR}/blood_vessel-edges -o ${DATA_DIR}/blood_vessel-csr

echo "[Convert cochlea graph]"
${REMOVE_SINGLETONS} -wf -i ${DATA_DIR}/cochlea-edges-uncompact -o ${DATA_DIR}/cochlea-edges
rm ${DATA_DIR}/cochlea-edges-uncompact
${EDGE_TO_CSR} -wf -i ${DATA_DIR}/cochlea-edges -o ${DATA_DIR}/cochlea-csr
