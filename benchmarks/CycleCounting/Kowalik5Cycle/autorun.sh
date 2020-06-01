#!/bin/bash
helpFunction()
{
   echo ""
   echo "Usage: $0 -g graph -r rounds -d orientation"
   echo -e "\t-g Path to graph input"
   echo -e "\t-r Number of rounds for each graph x param combination"
   echo -e "\t-d Type of orientation/edge directing"
   #echo -e "\t-c Description of what is parameterC"
   exit 1 # Exit script after printing help
}

while getopts "g:r:d:" opt
do
   case "$opt" in
      g ) graph="$OPTARG" ;;
      r ) rounds="$OPTARG" ;;
      d ) directType="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$graph" ] # || [ -z "$parameterB" ] || [ -z "$parameterC" ]
then
   echo "Graph parameter required";
   helpFunction
fi

if [ -z "$rounds" ]
then
    rounds="3"
fi

if [ -z "$directType" ]
then
    directType="0"
fi

# Begin script in case all parameters are correct
echo "$graph"
echo "$rounds"
echo "$directType"

# run the pure serial code
bazel run FiveCycle_main -- --serial -rounds "$rounds" -o "$directType" -s "$graph"
# run parallel with NUM_THREADS=1
export NUM_THREADS=1
bazel run FiveCycle_main -- -rounds "$rounds" -o "$directType" -s "$graph"
# run parallel with NUM_THREADS=72
export NUM_THREADS=72
bazel run FiveCycle_main -- -rounds "$rounds" -o "$directType" -s "$graph"
