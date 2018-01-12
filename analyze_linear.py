#!/usr/bin/python
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import argparse

def measure_error(filename):
  with open(filename,"r") as f:
    numeric = []
    for line in f:
      numeric.append(float(line))
    x = np.arange(float(len(numeric)))/float(len(numeric))
    analytic = np.zeros(len(numeric))
    analytic[np.where(x<0.75)] = 1.0

    L2 = np.mean((numeric - analytic)**2)
  return L2, x, numeric, analytic

if __name__== "__main__":
  parser = argparse.ArgumentParser(description='Plot convergence for fixed grid linear advection.')
  parser.add_argument('text_files',nargs='*')

  args = parser.parse_args()
  NX = []
  Error = []

  plt.figure()
  plt.ylabel("phi")
  plt.xlabel("x")
  plt.title("Lax-Friedrichs linear advection")
  x = np.arange(200)*0.005
  analytic = np.zeros(200)
  analytic[np.where(x<0.75)] = 1.0
  plt.plot(x,analytic,label='analytic')

  for filename in args.text_files:
    print filename
    resolution = filename.split('.')
    NX.append(float(resolution[1]))
    L2, x, numeric, analytic = measure_error(filename)
    Error.append(L2)
    plt.plot(x,numeric,'--',label='NX='+str(NX[-1]))


  plt.legend(loc='best')

  plt.figure()
  plt.title("Convergence plot")
  plt.ylabel("L_2 error")
  plt.xlabel("NX")
  plt.yscale('log')
  plt.xscale('log')
  plt.plot(NX,Error,'o')

  fit = np.polyfit(np.log(NX),np.log(Error),1)
  p = fit[0]
  A = np.exp(fit[1])
  plt.plot(NX, A * NX**p,label="NX^"+str(p))
  plt.legend(loc='best')
  plt.show()

