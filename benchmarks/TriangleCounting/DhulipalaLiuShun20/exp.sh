#!/bin/bash

let "e = $1 + 1000"
./Triangle -s -eg -blocksize 100 -n 2048 -bo $1 -be $e ../../../inputs/empty ../../../inputs/rMatGraph_J_2K_10K

# ./exp.sh 0 70
# ./exp.sh 71 140
