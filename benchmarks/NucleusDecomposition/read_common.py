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
      elif split[0].startswith("Data Structure Size"):
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