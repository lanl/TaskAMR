#!/usr/bin/python
import numpy as np
import subprocess
import sys
from analyze_linear import measure_error

def set_refinement_level(refinement_level):
  with open("global_const.rg","w") as f:
    f.write("-- required global constants\n")
    f.write("CELLS_PER_BLOCK_X = 2 -- must be multiple of 2\n")
    f.write("LEVEL_1_BLOCKS_X = 5\n")
    f.write("MAX_REFINEMENT_LEVEL = "+str(refinement_level)+"\n")
    f.write("NUM_PARTITIONS = 7\n")
    f.write("T_FINAL = 0.25\n")
    f.write("LENGTH_X = 1.0\n")
    f.close()

def test_single_resolution(refinement_level, filename, threshold, descriptor):
  ERROR = 0
  with open("/dev/null","w") as dev_null:

    set_refinement_level(refinement_level)
    subprocess.check_call(["../../github/legion/language/regent.py","1d_fix.rg"], stdout=dev_null)
    L2, x, numeric, analytic = measure_error(filename)
    if (L2 > threshold) or np.isnan(L2) :
      print descriptor+": \033[0;31mFAIL\033[0m ",L2," > ",threshold
      ERROR = 1
    else:
      print descriptor+": \033[0;32mPASS\033[0m ",L2," < ",threshold
  return ERROR

if __name__== "__main__":

  subprocess.check_call(["ln","-sf","linear_advection.rg","model.rg"])

  sys.exit(test_single_resolution(4, "linear.80.txt", 0.0487396, "Linear fixed NX=80") + test_single_resolution(5, "linear.160.txt", 0.0259696, "Linear fixed NX=160"))

