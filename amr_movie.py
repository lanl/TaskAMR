#!/usr/bin/python
#
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
# make movie with
# convert -delay 30 '*.png' movie.mov
# Or ImageJ File, Import, Image Sequence, File, SaveAs, AVI
#
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os

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

def save_fig(x, phi, time):

  fig, ax = plt.subplots()
  ax.plot(x, phi, '.', markersize=4)
  plt.ylim([-0.1,1.1])
  plt.savefig(time+'.png')

if __name__== "__main__":

  font = {'weight' : 'bold',
          'size' : 18}

  matplotlib.rc('font',**font)

  files = os.listdir("./")

  old_count = ''
  amr_files = []

  for file in files:
    texts = file.split("txt")
    if len(texts) > 1:
      entries = texts[0].split('_')
      print entries
      if entries[0] != old_count:
        if len(amr_files) > 0:
          x, phi = read_amr(amr_files)
          save_fig(x, phi, entries[0])
        old_count = entries[0]
        amr_files = []
        amr_files.append(file)
      else:
        amr_files.append(file)

#  plt.figure()
#  #plt.ylabel("phi")
#  #plt.xlabel("x")
#  #plt.title("AMR Lax-Friedrichs linear advection")
#  plt.ylim([-0.1,1.1])
#  plt.plot(x, phi, '.', markersize=8)
#  #plt.legend(loc='best')
  #plt.show()

