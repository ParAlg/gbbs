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
  program_pres = ["youtube", "orkut", "stackoverflow", "wiki", "ctr", "usa", "brain"]#["stackoverflow", "wiki", "ctr", "usa"]
  #empty = "empty_h"
  pred = ["ekcore"]#, "kcore"]
  num_rounds = 3
  e = 0.4
  d = 3
  batch_sizes = [1000000]#[10000, 100000, 1000000, 10000000]#[100, 1000, 10000,100000, 1000000, 10000000]
  num_workers = [60]#[1, 2, 4, 8, 16, 32, 60]
  #read_dir = "/home/jeshi/dynamic_graph/"
  write_dir = "/home/qliu19/exp-4-delete-dedup/"
  actual_batch_size = 1000000
  #for b_idx, b in enumerate(batch_sizes):
  for p_idx, p in enumerate(program_pres):
    for b_idx, b in enumerate(pred):
        read_filename = write_dir + str(b) + "_" + str(p) + "_" + str(e) + "_" + str(d) + "_" + str(1000000) + "_60.out"
        #num_lines = sum(1 for line in open(read_dir + pres + "_edges"))
        # if you are dynamic...round changes when you see ### Application...batch changes when you see batch running time
        # if you are static...
        best_avg = 0
        max_time = 0
        best_total = 0
        num_batches = 0
        with open(read_filename, "r") as read_file:
            cur_best = 0
            cur_max = 0
            for line in read_file:
                split = [x.strip() for x in line.split(':')]
                if "Application" in split[0]:
                    if best_total > cur_best or best_total == 0:
                        best_total = cur_best
                    if max_time < cur_max:
                        max_time = cur_max
                    cur_best = 0
                    cur_max = 0
                    num_batches = 0
                if split[0].startswith("### Batch Running Time"):
                    v = float(split[1])
                    if v > cur_max:
                        cur_max = v
                    cur_best += v
                    num_batches += 1
        best_avg = best_total/num_batches
        print(str(b), end = ",")
        print(str(p), end = ",")
        print(str(best_total), end = ",")
        print(str(best_avg), end = ",")
        print(str(max_time), end = "\n")

if __name__ == "__main__":
  main()
