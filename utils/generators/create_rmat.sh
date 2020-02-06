#!/bin/bash

#nums=(1048576 4194304 16777216 268435456 1073741824)
nums=(1048576)
edge_nums=(1 10 100 1000 10000 100000 1000000 2000000)
output_dir="/ssd0/graphs/edge_streams_coo"
for num in ${nums[@]}
do
  for edge_num in ${edge_nums[@]}
  do
    echo "numactl -i all ./RMAT -n $num -m $edge_num -outfile $output_dir/rmat_$num\_$edge_num"
    numactl -i all ./RMAT -n $num -m $edge_num -outfile $output_dir/rmat_$num\_$edge_num
  done
done
