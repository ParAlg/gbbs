import os
import sys
import signal
import time
import subprocess
from math import sqrt

def print_arr(row_headers, all_arr):
  for r_idx,row_head in enumerate(row_headers):
    arr = all_arr[r_idx]
    print(row_head + "\t"),
    for x in arr:
      print(str(x) + "\t"),
    print()

# "### Rank Running Time" Rank running time -- avg
  # "### Filter Graph Running Time" Filter running time -- avg
  # "### Table Running Time" Table + "### Count Running Time" Count + "### Peel Running Time" Peel running times -- avg and std dev
  # "rho" and "clique core"
  # "Data Structure Size"
# num_r_idx = number of r_idxs to process; num entries in file should be
# num_rounds * num_r_idx; append to end basically
def read_and_print(read_filename, num_rounds, avg_rank, avg_filter, avg_table,
  avg_count, avg_peel, std_dev_total, ds_size, r_idx, num_r_idx):
  rho = ""
  clique_core = ""
  data_struct_size = ""
  with open(read_filename, "r") as read_file:
    ranks = [[0]*num_rounds for x in range(num_r_idx)]
    filters = [[0]*num_rounds for x in range(num_r_idx)]
    tables = [[0]*num_rounds for x in range(num_r_idx)]
    counts = [[0]*num_rounds for x in range(num_r_idx)]
    peels = [[0]*num_rounds for x in range(num_r_idx)]
    idx = num_rounds - 1
    add_r_idx = 0
    first_round_flag = True
    for line in read_file:
      line = line.strip()
      split = [x.strip() for x in line.split(':')]
      if split[0].startswith("### Application"):
        idx += 1
        idx = idx % num_rounds
        if idx == 0 and not first_round_flag:
          add_r_idx += 1
        elif first_round_flag:
          first_round_flag = False
      elif split[0].startswith("### Rank Running Time"):
        ranks[add_r_idx][idx] = float(split[1])
      elif split[0].startswith("### Filter Graph Running Time"):
        filters[add_r_idx][idx] = float(split[1])
      elif split[0].startswith("### Table Running Time"):
        tables[add_r_idx][idx] = float(split[1])
      elif split[0].startswith("### Count Running Time"):
        counts[add_r_idx][idx] = float(split[1])
      elif split[0].startswith("### Peel Running Time"):
        peels[add_r_idx][idx] = float(split[1])
      elif split[0].startswith("rho"):
        rho = split[1]
      elif split[0].startswith("clique core"):
        clique_core = split[1]
      elif split[0].startswith("Data Structure Size") and ds_size[r_idx + add_r_idx] == "":
        ds_size[r_idx + add_r_idx] = split[1]
  for r in range(num_r_idx):
    t_ranks = ranks[r]
    t_filters = filters[r]
    t_tables = tables[r]
    t_counts = counts[r]
    t_peels = peels[r]
    avg_rank[r + r_idx] = sum(t_ranks) / num_rounds
    avg_filter[r + r_idx] = sum(t_filters) / num_rounds
    avg_table[r + r_idx] = sum(t_tables) / num_rounds
    avg_count[r + r_idx] = sum(t_counts) / num_rounds
    avg_peel[r + r_idx] = sum(t_peels) / num_rounds
    total_avg = avg_rank[r + r_idx] + avg_filter[r + r_idx] + avg_table[r + r_idx] + avg_count[r + r_idx] + avg_peel[r + r_idx]
    totals = [t_ranks[x] + t_filters[x] + t_tables[x] + t_counts[x] + t_peels[x] for x in range(num_rounds)]
    std_dev_total[r + r_idx] = sqrt(sum([(x - total_avg) ** 2 for x in totals]) / num_rounds)
  return (rho, clique_core)

def read_exp_1():
  r = "2"
  s = "3"
  pres = ["dblp", "as_skitter", "amazon", "youtube", "lj", "orkut", "friendster"]
  tts = ["1", "2", "3", "2", "3", "4", "5"]
  contigs = ["", "", "", "-contig", "-contig", "", ""]
  contig_pres = ["nc" if x == "" else "c" for x in contigs]
  relabel_pres = ["nr", "r"]
  efficients = ["0", "1", "2"]
  num_workers = [60]
  write_dir = "/home/jeshi_google_com/23_par_nd/"
  # We want one axis -- horiz -- to be over different tts, contigs, and nls
  # We want the other axis -- vert -- to be over different relabels and efficients
  # Print a new table for each file
  num_rounds = 4
  print("r: " + r + "s: " + s)
  for file_idx, filename in enumerate(pres):
    print("File: " + filename)
    rho = ""
    clique_core = ""
    data_struct_size = ""
    all_avg_rank = []
    all_avg_filter = []
    all_avg_table = []
    all_avg_count = []
    all_avg_peel = []
    all_std_dev_total = []
    all_ds_size = []
    row_headers = []
    for nw in num_workers:
      for relabel_pre in relabel_pres:
        for eff in efficients:
          # Here will begin a new line of the table
          avg_rank = []
          avg_filter = []
          avg_table = []
          avg_count = []
          avg_peel = []
          std_dev_total = []
          ds_size = []
          r_idx = 0
          for tt_idx, tt in enumerate(tts):
            contig = contigs[tt_idx]
            contig_pre = contig_pres[tt_idx]
            nls = ["2"]
            if (tt == "3" or tt == "4"):
              nls = [str(x) for x in range(2, int(r) + 1)]
            for nl in nls:
              avg_rank.append(0)
              avg_filter.append(0)
              avg_table.append(0)
              avg_count.append(0)
              avg_peel.append(0)
              std_dev_total.append(0)
              ds_size.append("")
            read_filename = write_dir + r + s + "_" + pres[file_idx] + "_" + tt + "_" + contig_pre + "_" + relabel_pre + "_" + eff +"_" + str(nw) + ".out"
            tup = read_and_print(read_filename, num_rounds, avg_rank, avg_filter, avg_table, avg_count, avg_peel, std_dev_total, ds_size, r_idx, len(nls))
            r_idx += len(nls)
            if (tup[0] != ""):
              rho = tup[0]
            if (tup[1] != ""):
              clique_core = tup[1]
          # now we must output
          all_avg_rank.append(avg_rank)
          all_avg_filter.append(avg_filter)
          all_avg_table.append(avg_table)
          all_avg_count.append(avg_count)
          all_avg_peel.append(avg_peel)
          all_std_dev_total.append(std_dev_total)
          all_ds_size.append(ds_size)
          row_headers.append(relabel_pre+eff)
    print("rho: " + rho + ", clique core: " + clique_core)
    print("Rank")
    print_arr(row_headers, all_avg_rank)
    print("\n Filter")
    print_arr(row_headers, all_avg_filter)
    print("\n Table")
    print_arr(row_headers, all_avg_table)
    print("\n Count")
    print_arr(row_headers, all_avg_count)
    print("\n Peel")
    print_arr(row_headers, all_avg_peel)
    print("\n Std Dev")
    print_arr(row_headers, all_std_dev_total)
    print("\n T Size")
    print_arr(row_headers, all_ds_size)
    print("\n\n")

def main():
  read_exp_1()

if __name__ == "__main__":
  main()