import os
import sys
import signal
import time
import subprocess

def signal_handler(signal,frame):
  print "bye\n"
  sys.exit(0)
signal.signal(signal.SIGINT,signal_handler)

def shellGetOutput(str) :
  process = subprocess.Popen(str,shell=True,stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
  output, err = process.communicate()


  #if (len(err) > 0):
  #    raise NameError(str+"\n"+output+err)
  return output

def computeTimeout(out):
  time = 0
  for line in out.splitlines():
    line = line.strip()
    split = [x.strip() for x in line.split(':')]
    if split[0].startswith("### Batch Running Time"):
      time += float(split[1])
  return time

def appendToFile(out, filename):
  with open(filename, "a+") as out_file:
    out_file.writelines(out)

# cluster time, total time, modularity, precision, recall, # comm,
# min comm, max comm, avg comm
def main():
  # Configured for Test 1
  program_dir = "../benchmarks/"
  programs = ["EdgeOrientation/ParallelLDS/LDS"] #["EdgeOrientation/LDS/LDS", "KCore/ApproximateKCore/KCore", "KCore/JulienneDBS17/KCore"]
  program_pres = ["plds"] #["lds", "kcore", "ekcore"]
  is_dynamic = [True]
  files = ["livejournal_deletion_edges"] #["youtube_deletion_edges","orkut_deletion_edges"]
  init_files = ["livejournal_insertion_edges"] #["youtube_insertion_edges", "orkut_insertion_edges"]#["livejournal_insertion_edges", "orkut_insertion_edges"]#["dblp_insertion_edges", "livejournal_insertion_edges"]
  pres = ["livejournal"] #["youtube_new", "orkut_new"]#["livejournal_1", "orkut_1"]#["dblp","livejournal"]
  empty = "empty_h"
  stats = "-stats"
  epss = [0.2, 0.4, 0.8, 1.6, 3.2, 6.4]
  deltas = [3, 6, 12, 24, 48, 96]
  batch_sizes = [1000000]#, 1000, 10000, 100000, 1000000, 10000000]
  num_workers = [60]#[1, 2, 4, 8, 16, 30, 60]
  read_dir = "/home/qliu19/dynamic_graph/"
  write_dir = "/home/qliu19/exp-1-delete/"
  for file_idx, filename in enumerate(files):
    num_lines = sum(1 for line in open(read_dir + filename))
    for program_idx, program in enumerate(programs):
      for e in epss:
        for d in deltas:
          for b in batch_sizes:
            for nw in num_workers:
              num_rounds = 3
              time = 0
              out_filename = write_dir + program_pres[program_idx] + "_" + pres[file_idx] + "_" + str(e) + "_" + str(d) + "_" + str(b) + "_" + str(nw) + ".out"
              batch_commands = []
              if is_dynamic[program_idx]:
                batch_commands = ["-b " + str(b)]
              elif b <= 10000:
                reasonable_b = 1000000
                batch_commands = batch_commands = ["-b "+str(b) +" -end_size "+str(x + reasonable_b) +" -start_size " + str(x) for x in range(0, num_lines, reasonable_b)]
              else:
                batch_commands = ["-start_size "+str(x)+" -num_dynamic_edges "+str(x + b) for x in range(0, num_lines, b)]
              for bc in batch_commands:
                ss = ("PARLAY_NUM_THREADS=" + str(nw) + " timeout 6h " + program_dir + program + " -s -i"
                " " + read_dir + filename + " -eps " + str(e) + " "
                "-delta " + str(d) + " " + bc + " "
                "-rounds " + str(num_rounds) + " -init_graph_file " + read_dir + init_files[file_idx]  + " " + stats + " " + read_dir + empty)
                print ss
                out = shellGetOutput(ss)
                appendToFile(out, out_filename)
                time += computeTimeout(out)
                if (time > 21600):
                  print("Timeout: " + ss)
                  break

if __name__ == "__main__":
  main()
