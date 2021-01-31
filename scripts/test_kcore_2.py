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

def appendToFile(out, filename):
  with open(filename, "a+") as out_file:
    out_file.writelines(out)

# cluster time, total time, modularity, precision, recall, # comm,
# min comm, max comm, avg comm
def main():
  # Configured for Test 1
  program_dir = "../benchmarks/"
  programs = ["EdgeOrientation/LDS/ParallelLDS"]
  #["EdgeOrientation/LDS/LDS", "KCore/ApproximateKCore/KCore", "KCore/JulienneDBS17/KCore"]
  program_pres = ["lds", "kcore", "ekcore"]
  is_dynamic = [True, False, False]
  files = ["dblp_edges","livejournal_edges"]
  pres = ["dblp","livejournal"]
  empty = "empty_h"
  stats = ""
  epss = [0.2, 0.4, 0.8, 1.6, 3.2, 6.4]
  deltas = [3, 6, 12, 24, 48, 96]
  batch_sizes = [100, 1000, 10000, 100000, 1000000, 10000000]
  num_workers = [60]#[1, 2, 4, 8, 16, 32, 60]
  read_dir = "/home/qliu19/dynamic_graph/"
  write_dir = "/home/qliu19/k-decomp-out/"
  for file_idx, filename in enumerate(files):
    for program_idx, program in enumerate(programs):
      for e in epss:
        for d in deltas:
          for b in batch_sizes:
            for nw in num_workers:
              out_filename = write_dir + program_pres[program_idx] + "_" + pres[file_idx] + "_" + str(e) + "_" + str(d) + "_" + str(b) + "_" + str(nw) + ".out"
              batch_commands = []
              if is_dynamic[program_idx]:
                batch_commands = ["-b " + str(b)]
              else:
                num_lines = sum(1 for line in open(read_dir + filename))
                batch_commands = ["-num_dynamic_edges " + str(x) for x in range(b, num_lines + 1, b)]
                if (num_lines % b != 0):
                  batch_commands.append(num_lines)
              for bc in batch_commands:
                ss = ("PARLAY_NUM_THREADS=" + str(nw) + " " + program_dir + program + " -s -i"
                " " + read_dir + filename + " -eps " + str(e) + " "
                "-delta " + str(d) + " " + bc + " "
                "-rounds 1 " + stats + " " + read_dir + empty)
                out = shellGetOutput(ss)
                appendToFile(out, out_filename)

if __name__ == "__main__":
  main()
