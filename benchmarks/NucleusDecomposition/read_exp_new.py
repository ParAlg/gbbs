import os
import sys
import signal
import time
import subprocess
from math import sqrt

from read_common import read_and_print, print_arr

def read_exp_new():
  r = "2"
  s = "3"
  pres = ["dblp", "as_skitter", "amazon", "youtube", "lj", "orkut", "friendster"]
  tts = ["2", "3", "4", "5"]
  contigs = ["-contig", "-contig", "", ""]
  contig_pres = ["nc" if x == "" else "c" for x in contigs]
  relabel_pres = ["r"]
  efficients = [ "1", "2"]
  num_workers = [60]
  write_dir = "/home/jeshi_google_com/23_par_nd/"
  # We want one axis -- horiz -- to be over different tts, contigs, and nls
  # We want the other axis -- vert -- to be over different relabels and efficients
  # Print a new table for each file
  num_rounds = 1
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
              nls = [str(x) for x in range(3, int(r) + 1)]
            for nl in nls:
              avg_rank.append(0)
              avg_filter.append(0)
              avg_table.append(0)
              avg_count.append(0)
              avg_peel.append(0)
              std_dev_total.append(0)
              ds_size.append("")
              read_filename = write_dir + r + s + "_" + pres[file_idx] + "_" + nl + "_" + tt + "_" + contig_pre + "_" + relabel_pre + "_" + eff +"_" + str(nw) + ".out"
              tup = read_and_print(read_filename, num_rounds, avg_rank, avg_filter, avg_table, avg_count, avg_peel, std_dev_total, ds_size, r_idx, 1)
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