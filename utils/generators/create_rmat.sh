#!/bin/bash


# nums=(32768)
# edge_nums=(32768000)
# nums=(2097152)
# edge_nums=(209715200)
# nums=(1048576)
# edge_nums=(1048576000)
nums=(16777216)
edge_nums=(1677721600)
output_dir="/home/sy/gbbs/utils/generators"
for num in ${nums[@]}
do
  for edge_num in ${edge_nums[@]}
  do
    echo "bazel run utils/generators:RMAT -- -n $num -m $edge_num -outfile $output_dir/rmat_$num\_$edge_num"
    bazel run utils/generators:RMAT -- -n $num -m $edge_num -outfile $output_dir/rmat_$num\_$edge_num
  done
done
