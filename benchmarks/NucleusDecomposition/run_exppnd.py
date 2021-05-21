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

def run_pnd():
  pres = ["dblp", "as_skitter", "amazon", "youtube", "lj", "orkut", "friendster"]
  algs = ["241", "240", "2400", "341", "340", "3400", "723", "734"]
  read_dir = "/home/jeshi_google_com/snap/"
  write_dir = "/home/jeshi_google_com/nd_exppnd/"
  for file_idx, filename in enumerate(pres):
    for alg in algs:
      for i in range(4):
        out_filename = write_dir + "pnd_" + filename + "_" + alg + ".out"
        ss = "./pnd " + read_dir + "out." + filename + " " + alg
        out = shellGetOutput(ss)
        appendToFile(out, out_filename)