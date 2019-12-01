#!/bin/bash
# experiments that use a base starting graph.

graphs=("/ssd1/graphs/bench_experiments/soc-LiveJournal1_sym.adj" "/ssd1/graphs/bench_experiments/com-orkut.ungraph.adj" "/ssd0/graphs/bench_experiments/usa_road.adj" "/ssd1/graphs/bench_experiments/twitter_sym.adj" "/ssd1/graphs/bench_experiments/friendster_sym.adj")

c_graphs=("/ssd0/graphs/bench_experiments/clueweb_sym.bytepda" "/ssd0/graphs/bench_experiments/hyperlink2014_sym.bytepda" "/ssd1/graphs/bench_experiments/hyperlink2012_sym.bytepda" )

#update_pcts=(0.2 0.4 0.6 0.8 1)
update_pcts=(1)
insert_to_queries=(0.25 0.5 0.75 1.0)
batch_sizes=(100000 1000000 10000000 100000000 1000000000 10000000000)

binaries=("jayanti.starting" "liutarjan.starting" "unite.starting" "unite_early.starting" "unite_nd.starting" "unite_rem_cas.starting" "unite_rem_lock.starting")

rounds=1
output_dir="data"
mkdir -p ${output_dir}
for graph in ${graphs[@]}
do
  NAME=$(basename $graph .adj)
  sync; echo 1 | sudo tee /proc/sys/vm/drop_caches
  for update_pct in ${update_pcts[@]}
  do
    for ins_to_query in ${insert_to_queries[@]}
    do
      for batch_size in ${batch_sizes[@]}
      do
        for binary in ${binaries[@]}
        do
          echo "sudo numactl -i all ./mains/${binary} -r ${rounds} -update_pct ${update_pct} -insert_to_query ${ins_to_query} -batch_size $batch_size $graph >> ${output_dir}/${NAME}_${update_pct}_${ins_to_query}_${batch_size}.data"
          sudo numactl -i all ./mains/${binary} -s -m -r ${rounds} -update_pct ${update_pct} -insert_to_query ${ins_to_query} -batch_size ${batch_size} $graph >> ${output_dir}/${NAME}_${update_pct}_${ins_to_query}_${batch_size}.data
        done
      done
    done
  done
done

for graph in ${c_graphs[@]}
do
  NAME=$(basename $graph .adj)
  sync; echo 1 | sudo tee /proc/sys/vm/drop_caches
  for update_pct in ${update_pcts[@]}
  do
    for ins_to_query in ${insert_to_queries[@]}
    do
      for batch_size in ${batch_sizes[@]}
      do
        for binary in ${binaries[@]}
        do
          echo "sudo numactl -i all ./mains/${binary} -c -s -m -r ${rounds} -update_pct ${update_pct} -insert_to_query ${ins_to_query} -batch_size $batch_size $graph >> ${output_dir}/${NAME}_${update_pct}_${ins_to_query}_${batch_size}.data"
          sudo numactl -i all ./mains/${binary} -c -s -m -r ${rounds} -update_pct ${update_pct} -insert_to_query ${ins_to_query} -batch_size ${batch_size} $graph >> ${output_dir}/${NAME}_${update_pct}_${ins_to_query}_${batch_size}.data
        done
      done
    done
  done
done
