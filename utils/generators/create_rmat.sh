#!/bin/bash

nums=(1048576 4194304 16777216 268435456 1073741824)
output_dir="/ssd1/graphs/edge_streams_coo"
for num in ${nums[@]}
do
  echo "numactl -i all ./RMAT -n $num -outfile $output_dir/rmat_$num"
  numactl -i all ./BarabasiAlbert -n $num -outfile $output_dir/rmat_$num
done
