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
import argparse

def trapezoid(x,f):
  value = np.sum(0.5 * (x[1:] - x[0:-1]) * (f[0:-1] + f[1:]))
  value += x[0] * f[0] + (1-x[-1])*f[-1]
  return value

def measure_error(filename):
  with open(filename,"r") as f:
    numeric = []
    for line in f:
      numeric.append(float(line))
    x = np.arange(float(len(numeric)))/float(len(numeric))
    x += 0.5 * x[1]
    analytic = np.zeros(len(numeric))
    analytic[np.where(x<0.75)] = 1.0

    L2 = trapezoid(x,(numeric-analytic)**2)
  return L2, x, numeric, analytic

if __name__== "__main__":
  import matplotlib
  import matplotlib.pyplot as plt

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
    print(filename)
    resolution = filename.split('.')
    NX.append(float(resolution[1]))
    L2, x, numeric, analytic = measure_error(filename)
    Error.append(L2)
    plt.plot(x,numeric,'--',label='NX='+str(NX[-1]))


  DX = 1.0 / np.array(NX)
  plt.legend(loc='best')

  plt.figure()
  plt.title("Convergence plot")
  plt.ylabel("L^2 error")
  plt.xlabel("dx")
  plt.ylim([5.0e-3,0.125])
  plt.xlim([1.0e-3,0.125])
  plt.yscale('log')
  plt.xscale('log')
  plt.plot(DX,Error,marker=(4,1), linewidth=0, markersize=15, label="Fixed", color='g')

  #fit = np.polyfit(np.log(DX[3:]),np.log(Error[3:]),1)
  #p = fit[0]
  #A = np.exp(fit[1])
  #plt.plot(DX, A * DX**p,linewidth=2,label="DX^"+ "%3.1f" % p)

  AMR_Error = [0.112082464854, 0.114230104721, 0.071889344041, 0.0485424093691, 0.0258501559003, 0.0126692353123, 0.00567836496466]
  #plt.plot(DX,AMR_Error,'.', linewidth=0, markersize=15, label="AMR", color='r')

  plt.legend(loc='best')

  font = {'family' : 'normal',
          'weight' : 'bold',
          'size' : 18}

  matplotlib.rc('font',**font)
  plt.show()

