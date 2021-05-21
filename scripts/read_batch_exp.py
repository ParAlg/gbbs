import os
import sys
import signal
import time
import subprocess
import numpy as np

def main():
  # Read each file in Test 1; separate LDS and Approx KCore
  # store average of average error, max of average error, total batch times, and average / max batch time
  # this is stored per epsilon and delta pair
  # Configured for Test 1 (1 round)
  program_dir = "../benchmarks/"
  programs = [ "EdgeOrientation/ParallelLDS/LDS"]# , "KCore/JulienneDBS17/KCore"]
  is_dynamic = [True]
  files = ["livejournal_edges"] #["livejournal_edges"]
  program_pres = [ "plds"]
  pres = ["livejournal"]
  empty = "empty_h"
  num_rounds = 100000
  e = [0.2, 0.4, 0.8, 1.6, 3.2, 6.4]
  d = [3, 6, 12, 24, 48, 96]
  batch_sizes = [1000000] # 1000000]#[100, 1000, 10000, 100000, 1000000, 10000000]
  num_workers = [60] #[1, 2, 4, 8, 16, 30, 60]
  read_dir = "/home/qliu19/dynamic_graph/"
  write_dir = "/home/qliu19/dynamic-graph-out-no-stats/"
  actual_batch_size = 1000000
  for file_name in pres:
    for eps in e:
        for l in d:
            for batch_size in batch_sizes:
                for worker in num_workers:
                    read_filename = write_dir + str(program_pres[0]) + "_" + file_name + "_" + str(eps) + "_" + str(l) + "_" + str(batch_size) + "_" + str(worker) + ".out"
                    #print(read_filename)
                    with open(read_filename, "r") as read_file:
                        round_idx = 0
                        prev_running_time = [0]*num_rounds
                        core_estimates = [0]*num_rounds
                        max_core_estimates = [0]*num_rounds
                        for line in read_file:
                            line = line.strip()
                            split = [x.strip() for x in line.split(':')]
                            if split[0].startswith("### Batch Running Time"):
                                prev_running_time[round_idx] = float(split[1])
                                round_idx += 1
                            elif split[0].startswith("### Per Vertex Average Coreness Error"):
                                core_estimates[round_idx] = float(split[1])
                            elif split[0].startswith("### Per Vertex Max Coreness Error"):
                                max_core_estimates[round_idx] = float(split[1])

                        max_runtime = max(prev_running_time)
                        total_runtime = sum(prev_running_time)
                        a = np.array(prev_running_time)
                        avg_runtime = a[np.nonzero(a)].mean()
                        min_runtime = min(i for i in prev_running_time if i > 0)
                        b = np.array(core_estimates)
                        avg_core_estimate = b[np.nonzero(b)].mean()
                        max_core_error = max(max_core_estimates)

                        print(file_name, end = ",")
                        print(str(eps), end = ",")
                        print(str(l), end = ",")
                        #print(avg_core_estimate, end = ",")
                        #print(max_core_error, end = ",")
                        #print(worker, end = ",")
                        print(total_runtime, end=",")
                        print(avg_runtime, end = ",")
                        print(max_runtime, end = "\n")

if __name__ == "__main__":
  main()
