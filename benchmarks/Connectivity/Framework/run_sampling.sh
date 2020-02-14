#!/bin/bash

graphs=("/ssd1/graphs/bench_experiments/soc-LiveJournal1_sym.adj" "/ssd1/graphs/bench_experiments/com-orkut.ungraph.adj" "/ssd0/graphs/bench_experiments/usa_road.adj" "/ssd1/graphs/bench_experiments/twitter_sym.adj" "/ssd1/graphs/bench_experiments/friendster_sym.adj")
c_graphs=("/ssd0/graphs/bench_experiments/clueweb_sym.bytepda" "/ssd0/graphs/bench_experiments/hyperlink2014_sym.bytepda" "/ssd1/graphs/bench_experiments/hyperlink2012_sym.bytepda" )

binaries=("benchmark_sampling")
output_dir="sampling_data"
rounds=5

mkdir -p ${output_dir}
mkdir -p ${output_dir}/json
mkdir -p ${output_dir}/raw

make benchmark_sampling

for graph in ${graphs[@]}
do
  NAME=$(basename $graph .adj)
  sync; echo 1 | sudo tee /proc/sys/vm/drop_caches
  for binary in ${binaries[@]}
  do
    echo "sudo numactl -i all ./${binary} -s -m -r ${rounds} $graph > ${output_dir}/raw/${NAME}.data"
    sudo numactl -i all ./${binary} -s -m -r ${rounds} $graph > ${output_dir}/raw/${NAME}.data
    echo "cat ${output_dir}/raw/${NAME}.data | grep \"^[^#;]\" > ${output_dir}/json/${NAME}.json"
    cat ${output_dir}/raw/${NAME}.data | grep "^[^#;]" > ${output_dir}/json/${NAME}.json
  done
done

for graph in ${c_graphs[@]}
do
  NAME=$(basename $graph .adj)
  sync; echo 1 | sudo tee /proc/sys/vm/drop_caches
  for binary in ${binaries[@]}
  do
    echo "sudo numactl -i all ./${binary} -s -c -m -r ${rounds} $graph > ${output_dir}/raw/${NAME}.data"
    sudo numactl -i all ./${binary} -s -c -m -r ${rounds} $graph > ${output_dir}/raw/${NAME}.data
    echo "cat ${output_dir}/raw/${NAME}.data | grep \"^[^#;]\" > ${output_dir}/json/${NAME}.json"
    cat ${output_dir}/raw/${NAME}.data | grep "^[^#;]" > ${output_dir}/json/${NAME}.json
  done
done
