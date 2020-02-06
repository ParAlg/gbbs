#!/bin/bash

#nums=(1048576 4194304 16777216 268435456 1073741824)
nums=(1073741824)
output_dir="/ssd0/graphs/edge_streams_coo"
for num in ${nums[@]}
do
  echo "numactl -i all ./BarabasiAlbert -n $num -edges_per_vertex 10 -outfile $output_dir/barabasialbert_10_$num"
  numactl -i all ./BarabasiAlbert -n $num -edges_per_vertex 10 -outfile $output_dir/barabasialbert_10_$num
done
