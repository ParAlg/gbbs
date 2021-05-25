import os
import sys
import signal
import time
import subprocess

def signal_handler(signal,frame):
  print "bye\n"
  sys.exit(0)
signal.signal(signal.SIGINT,signal_handler)

def shellGetOutput(str1) :
  process = subprocess.Popen(str1,shell=True,stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
  output, err = process.communicate()
  
  if (len(err) > 0):
    print(str1+"\n"+output+err)
  return output

def appendToFile(out, filename):
  with open(filename, "a+") as out_file:
    out_file.writelines(out)

def run_2():
  files = ["dblp_h", "as_skitter_h", "amazon_h", "youtube_h","lj_h", "orkut_h", "friendster_h"]
  pres = ["dblp", "as_skitter", "amazon", "youtube", "lj", "orkut", "friendster"]
  ss = ["6"]
  nl = "2"
  tts = ["1", "2", "3", "2", "3", "4", "5"]
  contigs = ["", "", "", "-contig", "-contig", "", ""]
  contig_pres = ["nc" if x == "" else "c" for x in contigs]
  relabels = ["", "-relabel"]
  efficients = ["0","1","2"]
  num_workers = [60]
  read_dir = "/home/jeshi/snap/"
  write_dir = "/home/jeshi/nd_exp1/"
  for s in ss:
    for file_idx, filename in enumerate(files):
      for tt_idx, tt in enumerate(tts):
        contig = contigs[tt_idx]
        contig_pre = contig_pres[tt_idx]
        for relabel in relabels:
          relabel_pre = "nr" if relabel == "" else "r"
          for eff in efficients:
            eff_pre = eff
            for nw in num_workers:
              for i in range(4):
                out_filename = write_dir + r + s + "_" + pres[file_idx] + "_" + nl + "_" tt + "_" + contig_pre + "_" + relabel_pre + "_" + eff_pre +"_" + str(nw) + ".out"
                ss = ("NUM_THREADS="+str(nw)+" timeout 6h bazel run :NucleusDecomposition_main -- "
                "-s -rounds 1 -r 2 -ss " + s + " -tt " + tt + " -nl"
                " " + nl + " " + relabel + " -efficient " + eff + " " + contig + " " + read_dir  + filename)
                out = shellGetOutput(ss)
                appendToFile(out, out_filename)


def run_3():
  r = "3"
  files = ["dblp_h", "as_skitter_h", "amazon_h", "youtube_h","lj_h", "orkut_h", "friendster_h"]
  pres = ["dblp", "as_skitter", "amazon", "youtube", "lj", "orkut", "friendster"]
  ss = ["6"]
  tts = ["1", "2", "3", "2", "3", "4", "5"]
  contigs = ["", "", "", "-contig", "-contig", "", ""]
  contig_pres = ["nc" if x == "" else "c" for x in contigs]
  relabels = ["", "-relabel"]
  efficients = ["0", "1","2"]
  num_workers = [60]
  read_dir = "/home/jeshi/snap/"
  write_dir = "/home/jeshi/nd_exp1/"
  for s in ss:
    for file_idx, filename in enumerate(files):
      for tt_idx, tt in enumerate(tts):
        contig = contigs[tt_idx]
        contig_pre = contig_pres[tt_idx]
        nls = ["2"]
        if (tt == "3" or tt == "4"):
          nls = ["2", "3"]
        for nl in nls:
          for relabel in relabels:
            relabel_pre = "nr" if relabel == "" else "r"
            for eff in efficients:
              eff_pre = eff
              for nw in num_workers:
                for i in range(4):
                  out_filename = write_dir + r + s + "_" + pres[file_idx] + "_" + nl + "_" tt + "_" + contig_pre + "_" + relabel_pre + "_" + eff_pre +"_" + str(nw) + ".out"
                  ss = ("NUM_THREADS="+str(nw)+" timeout 6h bazel run :NucleusDecomposition_main -- "
                  "-s -rounds 1 -r " + r + " -ss " + s + " -tt " + tt + " -nl"
                  " " + nl + " " + relabel + " -efficient " + eff + " " + contig + " " + read_dir  + filename)
                  out = shellGetOutput(ss)
                  appendToFile(out, out_filename)

def run_4():
  r = "4"
  files = ["dblp_h", "as_skitter_h", "amazon_h", "youtube_h","lj_h", "orkut_h", "friendster_h"]
  pres = ["dblp", "as_skitter", "amazon", "youtube", "lj", "orkut", "friendster"]
  ss = ["6"]
  tts = ["1", "2", "3", "2", "3", "4", "5"]
  contigs = ["", "", "", "-contig", "-contig", "", ""]
  contig_pres = ["nc" if x == "" else "c" for x in contigs]
  relabels = ["", "-relabel"]
  efficients = ["0", "1","2"]
  num_workers = [60]
  read_dir = "/home/jeshi/snap/"
  write_dir = "/home/jeshi/nd_exp1/"
  for s in ss:
    for file_idx, filename in enumerate(files):
      for tt_idx, tt in enumerate(tts):
        contig = contigs[tt_idx]
        contig_pre = contig_pres[tt_idx]
        nls = ["2"]
        if (tt == "3" or tt == "4"):
          nls = ["2", "3", "4"]
        for nl in nls:
          for relabel in relabels:
            relabel_pre = "nr" if relabel == "" else "r"
            for eff in efficients:
              eff_pre = eff
              for nw in num_workers:
                for i in range(4):
                  out_filename = write_dir + r + s + "_" + pres[file_idx] + "_" + nl + "_" tt + "_" + contig_pre + "_" + relabel_pre + "_" + eff_pre +"_" + str(nw) + ".out"
                  ss = ("NUM_THREADS="+str(nw)+" timeout 6h bazel run :NucleusDecomposition_main -- "
                  "-s -rounds 1 -r " + r + " -ss " + s + " -tt " + tt + " -nl"
                  " " + nl + " " + relabel + " -efficient " + eff + " " + contig + " " + read_dir  + filename)
                  out = shellGetOutput(ss)
                  appendToFile(out, out_filename)

def run_5():
  r = "5"
  files = ["dblp_h", "as_skitter_h", "amazon_h", "youtube_h","lj_h", "orkut_h", "friendster_h"]
  pres = ["dblp", "as_skitter", "amazon", "youtube", "lj", "orkut", "friendster"]
  ss = ["6"]
  tts = ["1", "2", "3", "2", "3", "4", "5"]
  contigs = ["", "", "", "-contig", "-contig", "", ""]
  contig_pres = ["nc" if x == "" else "c" for x in contigs]
  relabels = ["", "-relabel"]
  efficients = ["0", "1","2"]
  num_workers = [60]
  read_dir = "/home/jeshi/snap/"
  write_dir = "/home/jeshi/nd_exp1/"
  for s in ss:
    for file_idx, filename in enumerate(files):
      for tt_idx, tt in enumerate(tts):
        contig = contigs[tt_idx]
        contig_pre = contig_pres[tt_idx]
        nls = ["2"]
        if (tt == "3" or tt == "4"):
          nls = ["2", "3", "4", "5"]
        for nl in nls:
          for relabel in relabels:
            relabel_pre = "nr" if relabel == "" else "r"
            for eff in efficients:
              eff_pre = eff
              for nw in num_workers:
                for i in range(4):
                  out_filename = write_dir + r + s + "_" + pres[file_idx] + "_" + nl + "_" tt + "_" + contig_pre + "_" + relabel_pre + "_" + eff_pre +"_" + str(nw) + ".out"
                  ss = ("NUM_THREADS="+str(nw)+" timeout 6h bazel run :NucleusDecomposition_main -- "
                  "-s -rounds 1 -r " + r + " -ss " + s + " -tt " + tt + " -nl"
                  " " + nl + " " + relabel + " -efficient " + eff + " " + contig + " " + read_dir  + filename)
                  out = shellGetOutput(ss)
                  appendToFile(out, out_filename)