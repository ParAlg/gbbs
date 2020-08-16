#!/bin/bash

# let "bs  = 1000"
# let "s = $1 * $bs"
# let "e = $s + $bs"

# -blocksize 100
./Triangle -s -nb 10 -n 2048 -bo 0 ../../../inputs/empty ../../../inputs/rMatGraph_J_2K_10K

# ./Triangle -s -w -sg -trict 3 -nb  1 -n 7 -bo 0 ../../../inputs/small/triangles4.txt ../../../inputs/small/update4.txt

# ./Triangle -s -w -shuffle -nb  1 -n 7 -bo 0 ../../../inputs/small/triangles4.txt ../../../inputs/small/empty

# ./Triangle -s -blocksize 5 -n 128 -bo $1 ../../../inputs/empty ../../../inputs/rMatEdgesShuffle


# ./exp.sh 0 70
# ./exp.sh 71 140
