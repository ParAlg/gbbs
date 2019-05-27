#!/bin/bash

# Undirected algs
declare -a undir_alg=("CC" "Coloring" "DensestSubgraph" "KCore" "LDD" "MIS" "PageRank" "Spanner" "SpanningForest" "Triangle")
declare -a undir_mutable_alg=("Biconnectivity" "MaximalMatching" "SetCover")
declare -a unwgh_sssp=("BC" "BFS")
declare -a wgh_sssp=("BellmanFord" "wBFS" "WidestPath")
declare -a undir_wgh_mutable_alg=("MST")
declare -a dir_alg=("SCC")

declare undir_graph="/ssd1/graphs/bench_experiments/com-orkut.ungraph.adj"
declare unwgh_graph=${undir_graph}
declare wgh_graph="/ssd1/graphs/bench_experiments/com-orkut.ungraph_W.adj"

declare bench_path="../benchmark/"
declare mmap=true

declare datadir="timings/"
declare graphdir="${datadir}com-orkut/"

declare src=10012
declare rounds=3
declare numactl ="numactl -i all"

mkdir -p ${datadir}
mkdir -p ${graphdir}

# make the bench if necessary
pushd ${bench_path}; make -j; popd;

#export NUM_THREADS=1
unset NUM_THREADS
# Run undirected algorithms
for alg in "${undir_alg[@]}"; do
  declare extraflags=""
  if [ "$mmap" == "true" ]
  then
    extraflags+="-m "
  fi
  if [ "$alg" == "LDD" ]
  then
    extraflags+="-permute " # uses a random permutation in LDD
  fi
  if [ "$alg" == "KCore" ]
  then
    extraflags+="-nb 16 " # sets #buckets
  fi
  echo "Running ${alg}:"
  echo "cmd = \"${numactl} $bench_path${alg} -s -rounds ${rounds} ${extraflags} ${undir_graph} > ${graphdir}${alg}.dat\""
  ${numactl} $bench_path${alg} -s -rounds ${rounds} ${extraflags} ${undir_graph} > ${graphdir}${alg}.dat
done

# Run undirected mutable algorithms
for alg in "${undir_mutable_alg[@]}"; do
  declare extraflags=""
  if [ "$mmap" == "true" ]
  then
    extraflags+="-m "
  fi
  echo "Running ${alg}:"
  echo "cmd = \"${numactl} $bench_path${alg} -s -rounds 1 ${extraflags} ${undir_graph} > ${graphdir}${alg}.dat\""
  ${numactl} $bench_path${alg} -s -rounds 1 ${extraflags} ${undir_graph} > ${graphdir}${alg}.dat
  for i in {2 .. $rounds}
  do
    ${numactl} $bench_path${alg} -s -rounds 1 ${extraflags} ${undir_graph} >> ${graphdir}${alg}.dat
  done
done

# Run unweighted sssp algorithms
for alg in "${unwgh_sssp[@]}"; do
  declare extraflags=""
  if [ "$mmap" == "true" ]
  then
    extraflags+="-m "
  fi
  echo "Running ${alg}:"
  echo "cmd = \"${numactl} $bench_path${alg} -s -rounds ${rounds} -src ${src} ${extraflags} ${unwgh_graph} > ${graphdir}${alg}.dat\""
  ${numactl} $bench_path${alg} -s -rounds ${rounds} -src ${src} ${extraflags} ${unwgh_graph} > ${graphdir}${alg}.dat
done

# Run weighted sssp algorithms
for alg in "${wgh_sssp[@]}"; do
  declare extraflags=""
  if [ "$mmap" == "true" ]
  then
    extraflags+="-m "
  fi
  echo "Running ${alg}:"
  echo "cmd = \"${numactl} $bench_path${alg} -s -w -rounds ${rounds} -src ${src} ${extraflags} ${wgh_graph} > ${graphdir}${alg}.dat\""
  ${numactl} $bench_path${alg} -s -w -rounds ${rounds} -src ${src} ${extraflags} ${wgh_graph} > ${graphdir}${alg}.dat
done

# Run weighted mutable algorithms
for alg in "${undir_wgh_mutable_alg[@]}"; do
  declare extraflags=""
  if [ "$mmap" == "true" ]
  then
    extraflags+="-m "
  fi
  echo "Running ${alg}:"
  echo "cmd = \"${numactl} $bench_path${alg} -s -w -rounds 1 -src ${src} ${extraflags} ${wgh_graph} > ${graphdir}${alg}.dat\""
  ${numactl} $bench_path${alg} -s -w -rounds 1 -src ${src} ${extraflags} ${wgh_graph} > ${graphdir}${alg}.dat
  for i in {2 .. $rounds}
  do
    ${numactl} $bench_path${alg} -s -w -rounds 1 -src ${src} ${extraflags} ${wgh_graph} >> ${graphdir}${alg}.dat
  done
done
