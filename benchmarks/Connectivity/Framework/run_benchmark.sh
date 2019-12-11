#!/bin/bash

graphs=("/ssd1/graphs/bench_experiments/soc-LiveJournal1_sym.adj" "/ssd1/graphs/bench_experiments/com-orkut.ungraph.adj" "/ssd0/graphs/bench_experiments/usa_road.adj" "/ssd1/graphs/bench_experiments/twitter_sym.adj" "/ssd1/graphs/bench_experiments/friendster_sym.adj")
c_graphs=("/ssd0/graphs/bench_experiments/clueweb_sym.bytepda" "/ssd0/graphs/bench_experiments/hyperlink2014_sym.bytepda" "/ssd1/graphs/bench_experiments/hyperlink2012_sym.bytepda" )

binaries=("bfscc" "gbbscc" "jayanti_bfs" "jayanti_kout" "jayanti_ldd" "jayanti_nosample" "label_propagation" "liutarjan_nosample" "liutarjan_bfs" "liutarjan_ldd" "liutarjan_kout" "shiloach_vishkin" "unite_bfs" "unite_early_bfs" "unite_early_kout" "unite_early_ldd" "unite_early_nosample" "unite_kout" "unite_ldd" "unite_nd_bfs" "unite_nd_kout" "unite_nd_ldd" "unite_nd_nosample" "unite_nosample" "unite_rem_cas_bfs" "unite_rem_cas_kout" "unite_rem_cas_ldd" "unite_rem_cas_nosample" "unite_rem_lock_bfs" "unite_rem_lock_kout" "unite_rem_lock_ldd" "unite_rem_lock_nosample")
binary_dir="mains"

output_dir="data"

rounds=5

mkdir -p ${output_dir}
for graph in ${graphs[@]}
do
  NAME=$(basename $graph .adj)
  sync; echo 1 | sudo tee /proc/sys/vm/drop_caches
  rm ${output_dir}/${NAME}.data
  for binary in ${binaries[@]}
  do
    echo "sudo numactl -i all ./${binary_dir}/${binary} -s -m -r ${rounds} $graph >> ${output_dir}/${NAME}.data"
    sudo numactl -i all ./${binary_dir}/${binary} -s -m -r ${rounds} $graph >> ${output_dir}/${NAME}.data
  done
done

for graph in ${c_graphs[@]}
do
  NAME=$(basename $graph .adj)
  echo "sudo numactl -i all ./ConnectivityBenchmark -c -s -m -r 5 $graph > ${output_dir}/${NAME}.data"
  sync; echo 1 | sudo tee /proc/sys/vm/drop_caches
  rm ${output_dir}/${NAME}.data
  for binary in ${binaries[@]}
  do
    echo "sudo numactl -i all ./${binary_dir}/${binary} -c -s -m -r ${rounds} $graph >> ${output_dir}/${NAME}.data"
    sudo numactl -i all ./${binary_dir}/${binary} -c -s -m -r ${rounds} $graph >> ${output_dir}/${NAME}.data
  done
done
