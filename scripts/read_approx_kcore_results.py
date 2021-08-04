import os
import sys
import signal
import time
import subprocess

def signal_handler(signal,frame):
  print("bye")
  sys.exit(0)
signal.signal(signal.SIGINT,signal_handler)

def shellGetOutput(str) :
  process = subprocess.Popen(str,shell=True,stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
  output, err = process.communicate()

  if (len(err) > 0):
    raise NameError(str+"\n"+output.decode('utf-8')+err.decode('utf-8'))
  return output.decode('utf-8')

def appendToFile(out, filename):
  with open(filename, "a+") as out_file:
    out_file.writelines(out)

def benchmarkToProgramPath(benchmark):
  benchmark_dict = {
    "LDS" : "EdgeOrientation/LDS/LDS",
    "ParallelLDS" : "EdgeOrientation/ParallelLDS/LDS",
    "KCore" : "KCore/ApproximateKCore/KCore"
  }
  return benchmark_dict.get(benchmark)

def benchmarkToIsDynamic(benchmark):
  benchmark_dict = {
    "LDS" : True,
    "ParallelLDS" : True,
    "KCore" : False
  }
  return benchmark_dict.get(benchmark)

def main():
  git_init_process = subprocess.Popen("git init ..",shell=True,stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
  output, err = git_init_process.communicate()

  # Read parameters from setup file
  with open('approx_kcore_setup.txt') as parameters_file:
    for line in parameters_file:
      line = line.strip()
      split = [x.strip() for x in line.split(':')]
      if len(split) <= 1:
        continue
      params = [x.strip() for x in split[1].split(',')]
      if line.startswith("Input graph directory"):
        read_dir = split[1]
      elif line.startswith("Dynamic graphs"):
        files = params.copy()
      elif line.startswith("Output directory"):
        write_dir = split[1]
      elif line.startswith("Benchmarks"):
        programs = [benchmarkToProgramPath(x) for x in params]
        is_dynamic = [benchmarkToIsDynamic(x) for x in params]
        program_pres = params.copy()
      elif line.startswith("Numbers of workers"):
        num_workers = params.copy()
      elif line.startswith("Epsilons"):
        epss = params.copy()
      elif line.startswith("Lambdas/Deltas"):
        deltas = params.copy()
      elif line.startswith("Number of Levels Divisor"):
        divisors = params.copy()
      elif line.startswith("Batch sizes"):
        batch_sizes = params.copy()
      elif line.startswith("Output stats"):
        if split[1] == "True":
          stats = "-stats"
        else:
          stats = ""
      elif line.startswith("Output sizes"):
        if split[1] == "True":
          size = "-size"
        else:
          size = ""
      elif line.startswith("Opt"):
        if split[1] == "True":
          opt = "-opt"
        else:
          opt = ""
  # Setup other parameters
  program_dir = "../benchmarks/"
  empty = "../benchmarks/EdgeOrientation/ParallelLDS/empty_h"
  for file_idx, filename in enumerate(files):
    for program_idx, program in enumerate(programs):
      for e in epss:
        for d in deltas:
          for divisor in divisors:
            for b in batch_sizes:
              for nw in num_workers:
                num_rounds = 4
                out_path_components = [program_pres[program_idx], filename, e,
                        d, b, nw, divisor, ".out"]
                read_filename = os.path.join(write_dir, "_".join(out_path_components))

                best_avg_time = 0
                best_max_time = 0
                best_total_time = 0
                best_avg_error = 0
                best_max_error = 0

                with open(read_filename, "r") as read_file:
                    cur_max_time = 0
                    cur_total_time = 0
                    num_iterations = 0
                    cur_max_error = 0
                    cur_total_error = 0
                    for line in read_file:
                        line = line.strip()
                        if "------------" in line and cur_max_time > 0:
                            cur_avg_time = cur_total_time / num_iterations
                            if best_avg_time == 0 or cur_avg_time < best_avg_time:
                                best_avg_time = cur_avg_time
                                best_max_time = cur_max_time
                                best_total_time = cur_total_time
                                best_avg_error = cur_total_error / num_iterations
                                best_max_error = cur_max_error
                            cur_max_time = 0
                            cur_total_time = 0
                            num_iterations = 0
                            cur_max_error = 0
                            cur_total_error = 0
                        else:
                            split = [x.strip() for x in line.split(':')]
                            if split[0].startswith("### Batch Running Time"):
                                cur_total_time += float(split[1])
                                if float(split[1]) > cur_max_time:
                                    cur_max_time = float(split[1])
                                num_iterations += 1
                            elif split[0].startswith("### Per Vertex Average"):
                                cur_total_error += float(split[1])
                            elif split[0].startswith("### Per Vertex Max"):
                                if float(split[1]) > cur_max_error:
                                    cur_max_error = float(split[1])

                    cur_avg_time = cur_total_time / num_iterations
                    if best_avg_time == 0 or cur_avg_time < best_avg_time:
                        best_avg_time = cur_avg_time
                        max_time = cur_max_time
                        best_total_time = cur_total_time
                        best_avg_error = cur_total_error / num_iterations
                        best_max_error = cur_max_error
                for param in out_path_components:
                    print(str(param), end = ",")
                print(str(best_avg_time), end = ",")
                print(str(best_max_time), end = ",")
                print(str(best_total_time), end = ",")
                print(str(best_avg_error), end=",")
                print(str(best_max_error), end="\n")

if __name__ == "__main__":
  main()
