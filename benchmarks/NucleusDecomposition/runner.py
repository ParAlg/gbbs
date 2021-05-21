import os
import sys
import signal
import time
import subprocess
from run_exp_1 import run_2, run_3, run_4, run_ktruss
from run_exppnd import run_pnd

def main():
  run_2()
  run_3()
  run_4()

if __name__ == "__main__":
  main()
