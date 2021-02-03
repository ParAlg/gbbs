import os
import sys
import signal
import time
import subprocess

def main():
  # Read each file in Test 1; separate LDS and Approx KCore
  # store average of average error, max of average error, total batch times, and average / max batch time
  # this is stored per epsilon and delta pair
  # Configured for Test 1 (1 round)
  #program_dir = "../benchmarks/"
  #programs = [ "KCore/ApproximateKCore/KCore"]# , "KCore/JulienneDBS17/KCore"]
  #is_dynamic = [True, False, False]
  #files = ["dblp_edges", "livejournal_edges"]
  program_pres = [ "lds"]
  pres = "dblp_insertion_perf"
  #empty = "empty_h"
  num_rounds = 2
  e = 6.4
  d = 3
  batch_sizes = [ 100000]#100, 1000, 10000,, 1000000, 10000000
  num_workers = [60]#[1, 2, 4, 8, 16, 32, 60]
  #read_dir = "/home/jeshi/dynamic_graph/"
  write_dir = "/home/sy/dogfood-out/"
  actual_batch_size = 100000
  total_batch_times = [0]*num_rounds
  avg_batch_times = [0]*num_rounds
  max_batch_times = [0]*num_rounds
  #for b_idx, b in enumerate(batch_sizes):
  for p_idx, p in enumerate(program_pres):
    print("", end = "\n")
    print(p, end = "\n")
    for b_idx, b in enumerate(batch_sizes):
      read_filename = write_dir + str(p) + "_" + pres + "_" + str(e) + "_" + str(d) + "_" + str(b) + "_60.out"
    #num_lines = sum(1 for line in open(read_dir + pres + "_edges"))
    # if you are dynamic...round changes when you see ### Application...batch changes when you see batch running time
    # if you are static...
      with open(read_filename, "r") as read_file:
        sum_runtime = [0]*num_rounds
        max_runtime = [0]*num_rounds
        num_batches = [0]*num_rounds
        round_idx = num_rounds - 1
        prev_running_time = [0]*num_rounds
        for line in read_file:
          line = line.strip()
          split = [x.strip() for x in line.split(':')]
          if split[0].startswith("### Application"): #and is_dynamic[p_idx]:#and split[1].startswith("LDS"):
            round_idx = round_idx + 1
            round_idx = round_idx % num_rounds
          if split[0].startswith("### Batch Running Time"):
            prev_running_time[round_idx] += float(split[1])
          elif split[0].startswith("### Batch Num"):
            #if int(split[1]) % b == 0 or int(split[1]) == num_lines:
            num_batches[round_idx] += 1
            sum_runtime[round_idx] += prev_running_time[round_idx]
            max_runtime[round_idx] = max(max_runtime[round_idx], prev_running_time[round_idx])
            prev_running_time[round_idx] = 0
        for i in range(num_rounds):
          total_batch_times[i] = sum_runtime[i]
          avg_batch_times[i] = sum_runtime[i] if num_batches[i] == 0 else sum_runtime[i] / num_batches[i]
          max_batch_times[i] = max_runtime[i]
      # now we must output
      min_idx = 1
      min_total_batch_times = total_batch_times[1]
      for i in range(1, num_rounds):
        if (total_batch_times[i] < min_total_batch_times):
          min_total_batch_times = total_batch_times[i]
          min_idx = i
      print(str(b), end = "\t")
      print(total_batch_times[min_idx], end = "\t")
      print(avg_batch_times[min_idx], end = "\t")
      print(max_batch_times[min_idx], end = "\n")
      for i in range(num_rounds):
        total_batch_times[i] = 0
        avg_batch_times[i] = 0
        max_batch_times[i] = 0

if __name__ == "__main__":
  main()
