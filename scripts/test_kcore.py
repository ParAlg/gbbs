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
  programs = ["EdgeOrientation/LDS/LDS", "KCore/ApproximateKCore/KCore", "KCore/JulienneDBS17/KCore"]
  program_pres = ["lds", "kcore", "ekcore"]
  is_dynamic = [True, False, False]
  files = ["dblp_edges","livejournal_edges"]
  pres = ["dblp","livejournal"]
  empty = "empty_h"
  stats = ""
  epss = [0.4] #[0.2, 0.4, 0.8, 1.6, 3.2, 6.4]
  deltas = [3] #[3, 6, 12, 24, 48, 96]
  batch_sizes = [100]#, 1000, 10000, 100000, 1000000, 10000000]
  num_workers = [60]#[1, 2, 4, 8, 16, 32, 60]
  read_dir = "/home/jeshi/dynamic_graph/"
  write_dir = "/home/jeshi/dogfood-out/"
  for file_idx, filename in enumerate(files):
    for program_idx, program in enumerate(programs):
      for e in epss:
        for d in deltas:
          for b in batch_sizes:
            for nw in num_workers:
              time = 0
              out_filename = write_dir + program_pres[program_idx] + "_" + pres[file_idx] + "_" + str(e) + "_" + str(d) + "_" + str(b) + "_" + str(nw) + ".out"
              batch_commands = []
              if is_dynamic[program_idx]:
                batch_commands = ["-b " + str(b)]
              else:
                num_lines = sum(1 for line in open(read_dir + filename))
                #batch_commands = ["-b "+str(b) +" -end_size "+str(x + 100*b) +" -start_size " + str(x) for x in range(0, num_lines, 100*b)]
                batch_commands = ["-start_size "+str(x)+" -num_dynamic_edges "+str(x + b) for x in range(0, num_lines, b)]
                #if (num_lines % b != 0):
                #  batch_commands.append("-num_dynamic_edges " + str(num_lines))
              for bc in batch_commands:
                ss = ("PARLAY_NUM_THREADS=" + str(nw) + " timeout 6h " + program_dir + program + " -s -i"
                " " + read_dir + filename + " -eps " + str(e) + " "
                "-delta " + str(d) + " " + bc + " "
                "-rounds 4 " + stats + " " + read_dir + empty)
                out = shellGetOutput(ss)
                appendToFile(out, out_filename)
                time += computeTimeout(out)
                if (time > 21600):
                  print("Timeout: " + ss)
                  break

if __name__ == "__main__":
  main()
