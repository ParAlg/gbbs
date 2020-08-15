#!/bin/bash

let "e = $1 + 100"
./Triangle -s -eg -blocksize 100 -n 256 -bo $1 -be $e ../../../inputs/empty ../../../inputs/rMatGraph_J_200_1K_300

# ./Triangle -s -eg -blocksize 5 -n 128 -bo $1 -be $e ../../../inputs/empty ../../../inputs/rMatEdgesShuffle


# ./exp.sh 0 70
# ./exp.sh 71 140
