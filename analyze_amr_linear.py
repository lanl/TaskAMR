#!/usr/bin/python
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import argparse

if __name__== "__main__":
  parser = argparse.ArgumentParser(description='Plot convergence for fixed grid linear advection.')
  parser.add_argument('text_files',nargs='*')

  args = parser.parse_args()
  x = []
  phi = []


  for filename in args.text_files:
    with open(filename, "r") as f:
      for line in f:
        data = line.split()
        x.append(data[0])
        phi.append(data[1])

plt.figure()
plt.ylabel("phi")
plt.xlabel("x")
plt.title("AMR Lax-Friedrichs linear advection")
plt.ylim([-0.1,1.1])
plt.plot(x, phi, '.', label='numeric')
plt.legend(loc='best')
plt.show()

