#!/usr/bin/env python3
#
# Copyright (c) 2018, Triad National Security, LLC
# All rights reserved.
# 
# This program was produced under U.S. Government contract 89233218CNA000001 for
# Los Alamos National Laboratory (LANL), which is operated by Triad National
# Security, LLC for the U.S. Department of Energy/National Nuclear Security
# Administration.
# 
# THIS SOFTWARE IS PROVIDED BY TRIAD NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL TRIAD NATIONAL SECURITY, LLC OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
# 
# If software is modified to produce derivative works, such modified software should be
# clearly marked, so as not to confuse it with the version available from LANL.
import numpy as np
import os
import subprocess
import sys
from analyze_linear import measure_error

legion_root = os.environ.get('LEGION_ROOT', '../../github/legion')
regent = os.path.join(legion_root, 'language/regent.py')

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
    subprocess.check_call([regent,"1d_fix.rg"], stdout=dev_null)
    L2, x, numeric, analytic = measure_error(filename)
    if (L2 > threshold) or np.isnan(L2) :
      print(descriptor+": \033[0;31mFAIL\033[0m ",L2," > ",threshold)
      ERROR = 1
    else:
      print(descriptor+": \033[0;32mPASS\033[0m ",L2," < ",threshold)
  return ERROR

if __name__== "__main__":

  subprocess.check_call(["ln","-sf","linear_advection.rg","model.rg"])

  sys.exit(test_single_resolution(4, "linear.80.txt", 0.0487396, "Linear fixed NX=80") + test_single_resolution(5, "linear.160.txt", 0.0259696, "Linear fixed NX=160"))

