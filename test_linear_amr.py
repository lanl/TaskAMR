#!/usr/bin/python
import numpy as np
import subprocess
import sys
from analyze_linear import trapezoid
from analyze_amr_linear import read_amr
from test_linear import set_refinement_level

def test_amr(refinement_level, filenames, threshold, descriptor):
  ERROR = 0
  with open("/dev/null","w") as dev_null:

    set_refinement_level(refinement_level)
    subprocess.check_call(["../../github/legion/language/regent.py","1d_amr.rg"], stdout=dev_null)

    x, numeric = read_amr(filenames)
    x = np.array(x)
    analytic = np.zeros(len(x))
    analytic[np.where(np.array(x)<0.75)] = 1.0
    L2 = trapezoid(x, (numeric - analytic)**2)
 
    if (L2 > threshold) or np.isnan(L2) :
      print descriptor+": \033[0;31mFAIL\033[0m ",L2," > ",threshold
      ERROR = 1
    else:
      print descriptor+": \033[0;32mPASS\033[0m",L2," < ",threshold
  return ERROR

if __name__== "__main__":

  subprocess.check_call(["ln","-sf","linear_advection.rg","model.rg"])
  subprocess.check_call(["ln","-sf","linear_advection_amr.rg","model_amr.rg"])

  sys.exit(test_amr(4, ["linear_amr.20.0.txt","linear_amr.40.0.txt",
           "linear_amr.80.0.txt",
           "linear_amr.80.14.txt","linear_amr.80.28.txt"], 0.0487396, "AMR 4 levels"))
