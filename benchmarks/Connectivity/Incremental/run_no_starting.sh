#!/bin/bash
# experiments that use a base starting graph.

graphs=("/ssd0/graphs/edge_streams_coo/barabasialbert_10_1048576.coo" "/ssd0/graphs/edge_streams_coo/barabasialbert_10_4194304.coo" "/ssd0/graphs/edge_streams_coo/barabasialbert_10_16777216.coo" "/ssd0/graphs/edge_streams_coo/barabasialbert_10_67108864.coo" "/ssd0/graphs/edge_streams_coo/barabasialbert_10_268435456.coo" "/ssd0/graphs/edge_streams_coo/barabasialbert_10_1073741824.coo" "/ssd1/graphs/edge_streams_coo/rmat_1048576" "/ssd1/graphs/edge_streams_coo/rmat_4194304" "/ssd1/graphs/edge_streams_coo/rmat_16777216" "/ssd1/graphs/edge_streams_coo/rmat_268435456")

update_pcts=(1)
insert_to_queries=(2)
batch_sizes=(1000 10000 100000 1000000 10000000 100000000 1000000000)

binaries=("jayanti.no_starting" "liutarjan.no_starting" "unite.no_starting" "unite_early.no_starting" "unite_nd.no_starting" "unite_rem_cas.no_starting" "unite_rem_lock.no_starting" )

rounds=1

output_dir="no_starting"
mkdir -p ${output_dir}
for graph in ${graphs[@]}
do
  NAME=$(basename $graph .coo)
    for ins_to_query in ${insert_to_queries[@]}
    do
      for batch_size in ${batch_sizes[@]}
      do
        for binary in ${binaries[@]}
        do
          echo "sudo numactl -i all ./mains/${binary} -r ${rounds} -insert_to_query ${ins_to_query} -batch_size $batch_size -in_file ${graph}  >> ${output_dir}/basic/${NAME}_${ins_to_query}_${batch_size}.data"
          sudo numactl -i all ./mains/${binary} -r ${rounds} -insert_to_query ${ins_to_query} -batch_size $batch_size -in_file ${graph}  >> ${output_dir}/${NAME}_${ins_to_query}_${batch_size}.data
        done
      done
    done
done
