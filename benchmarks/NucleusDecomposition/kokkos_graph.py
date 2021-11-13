import os
import sys
import signal
import time
import subprocess

def main():
  read_dir = "/home/jeshi/snap/"
  grpah_name = "com-dblp.ungraph.txt"
  graph_pre = "dblp"
  read_filename = read_dir + graph_name
  write_filename = read_dir + graph_pre + ".kokkos"
  with open(read_filename, "r") as read_file:
    with open(write_filename, "w") as write_file:
      for line in read_file:
        line = line.strip()
        split = [x.strip() for x in line.split('\t')]
        write_file.write(split[0] + "\t" + split[1] + "\n")
        write_file.write(split[1] + "\t" + split[0] + "\n")

if __name__ == "__main__":
  main()