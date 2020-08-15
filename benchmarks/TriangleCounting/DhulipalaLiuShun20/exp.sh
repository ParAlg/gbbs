#!/bin/bash

# let "bs  = 1000"
# let "s = $1 * $bs"
# let "e = $s + $bs"

# -blocksize 100
./Triangle -s -eg -nb 10 -n 2048 -bo $1 ../../../inputs/empty ../../../inputs/rMatGraph_J_2K_10K

# ./Triangle -s -eg -blocksize 5 -n 128 -bo $1 -be $e ../../../inputs/empty ../../../inputs/rMatEdgesShuffle


# ./exp.sh 0 70
# ./exp.sh 71 140
