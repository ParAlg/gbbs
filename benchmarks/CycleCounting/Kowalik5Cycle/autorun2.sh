helpFunction()
{
   echo ""
   echo "Usage: $0 -g graph -f file"
   echo -e "\t-g Path to graph input"
   echo -e "\t-f Path to output file"
   exit 1;
}

while getopts "g:f:" opt
do
   case "$opt" in
      g ) graph="$OPTARG" ;;
      f ) file="$OPTARG" ;;
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$graph" ] || [ -z "$file" ] # || [ -z "$parameterC" ]
then 
   echo "Graph and file parameter required";
   helpFunction
fi

echo "Running experiment on $graph, output to $file"


for ((directType = 0; directType < 3; directType++)) # 3X
do for ((n = 100; n <= 100000000; n = n<<2)) # 5X
    do
        export NUM_THREADS=1
        timeout 18000 bazel run FiveCycle_main -- --serial -rounds 3 -b "$n" -o "$directType" -s "$graph"
        for ((thr = 1; thr <= 72; thr = thr << 1)) # 7X
        do
        export NUM_THREADS=$thr
        timeout 18000 bazel run FiveCycle_main -- -rounds 3 -b "$n" -o "$directType" -s "$graph"
        timeout 18000 bazel run FiveCycle_main -- --no-schedule -rounds 3 -b "$n" -o "$directType" -s "$graph"
        done
        export NUM_THREADS=72
        timeout 18000 bazel run FiveCycle_main -- -rounds 3 -b "$n" -o "$directType" -s "$graph"
        timeout 18000 bazel run FiveCycle_main -- --no-schedule -rounds 3 -b "$n" -o "$directType" -s "$graph"
   done
done | python3 output-reader.py "$file"