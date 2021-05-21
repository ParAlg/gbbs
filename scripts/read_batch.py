import os
import sys
import signal
import time
import subprocess
import numpy as np

def main():
    program_dir = "../benchmarks/"
    programs = [ "EdgeOrientation/ParallelLDS/LDS"]# , "KCore/JulienneDBS17/KCore"]
    is_dynamic = [True]
    files = ["dblp_edges", "livejournal_edges"]
    program_pres = [ "plds"]
    pres = ["orkut-2"]
    empty = "empty_h"
    num_rounds = 10000000
    e = 0.4
    d = 3
    batch_sizes = [1000000]#[100, 1000, 10000, 100000, 1000000, 10000000]
    num_workers = [60]#[1, 2, 4, 8, 16, 32, 60]
    read_dir = "/home/qliu19/dynamic_graph/"
    write_dir = "/home/qliu19/dynamic-graph-out-large/"
    actual_batch_size = 1000000
    total_batch_times = [0]*num_rounds
    avg_batch_times = [0]*num_rounds
    max_batch_times = [0]*num_rounds
    for file_name in pres:
        for batch in batch_sizes:
            read_filename = write_dir + str(program_pres[0]) + "_" + file_name + "_" + str(e) + "_" + str(d) + "_" + str(batch) + "_60.out"
            with open(read_filename, "r") as read_file:
                round_idx = 0
                prev_running_time = [0]*num_rounds
                for line in read_file:
                    line = line.strip()
                    split = [x.strip() for x in line.split(':')]
                    if split[0].startswith("### Batch Running Time"):
                        prev_running_time[round_idx] = float(split[1])
                        round_idx += 1
            max_runtime = max(prev_running_time)
            a = np.array(prev_running_time)
            avg_runtime = a[np.nonzero(a)].mean()
            min_runtime = min(i for i in prev_running_time if i > 0)
            total_runtime = sum(a)/3

            print(file_name, end = ",")
            #print(batch, end=",")
            print(total_runtime, end=",")
            print(avg_runtime, end = ",")
            print(max_runtime, end = "\n")

if __name__ == "__main__":
    main()
