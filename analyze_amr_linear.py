#!/usr/bin/python
# Copyright (c) 2018, Los Alamos National Security, LLC
# All rights reserved.
#
# Copyright 2018. Los Alamos National Security, LLC. This software was produced under
# U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL),
# which is operated by Los Alamos National Security, LLC for the U.S. Department of
# Energy. The U.S. Government has rights to use, reproduce, and distribute this
# software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY
# WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
# If software is modified to produce derivative works, such modified software should be
# clearly marked, so as not to confuse it with the version available from LANL.
#
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import argparse

def read_amr(filenames):
  x = []
  phi = []

  for filename in filenames:
    with open(filename, "r") as f:
      for line in f:
        data = line.split()
        x.append(float(data[0]))
        phi.append(float(data[1]))

  return x,phi
  
if __name__== "__main__":
  parser = argparse.ArgumentParser(description='Plot convergence for fixed grid linear advection.')
  parser.add_argument('text_files',nargs='*')

  args = parser.parse_args()

  x, phi = read_amr(args.text_files)

  plt.figure()
  plt.ylabel("phi")
  plt.xlabel("x")
  plt.title("AMR Lax-Friedrichs linear advection")
  plt.ylim([-0.1,1.1])
  plt.xlim([-0.1,1.1])
  plt.plot(x, phi, '.', label='numeric')
  plt.legend(loc='best')
  plt.show()

