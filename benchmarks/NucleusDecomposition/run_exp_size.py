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


def run_size():
  files = ["dblp_h", "amazon_h", "youtube_h","as_skitter_h", "lj_h", "orkut_h", "friendster_h"]
  pres = ["dblp", "amazon", "youtube", "as_skitter", "lj", "orkut", "friendster"]
  sss = ["3", "4", "4", "5"]
  rs = ["2","2","3","4"]
  tts = ["3", "4"]
  contigs = ["-contig", ""]
  contig_pres = ["nc" if x == "" else "c" for x in contigs]
  relabels = ["-relabel"]
  efficients = ["1"]
  num_workers = [60]
  read_dir = "/home/jeshi_google_com/snap/"
  write_dir = "/home/jeshi_google_com/nd_exp1/"
  for file_idx, filename in reversed(list(enumerate(files))):
    for s_idx, s in enumerate(sss):
      r = rs[s_idx]
      if (s_idx == 3) and (file_idx == 6):
        continue
      for tt_idx, tt in enumerate(tts):
        contig = contigs[tt_idx]
        contig_pre = contig_pres[tt_idx]
        nls = [str(x) for x in range(2, int(r) + 1)]
        for nl in nls:
          for relabel in relabels:
            relabel_pre = "nr" if relabel == "" else "r"
            for eff in efficients:
              eff_pre = eff
              for nw in num_workers:
                out_filename = write_dir + r + s + "_" + pres[file_idx] + "_" + nl + "_" + tt + "_" + contig_pre + "_" + relabel_pre + "_" + eff_pre +"_" + str(nw) + ".out"
                ss = ("NUM_THREADS="+str(nw)+" timeout 6h bazel run :NucleusDecomposition_main -- "
                "-s -rounds 1 -r " + r + " -ss " + s + " -tt " + tt + " -output_size -nl "
                " " + nl + " " + relabel + " -efficient " + eff + " " + contig + " " + read_dir  + filename)
                out = shellGetOutput(ss)
                appendToFile(out, out_filename)