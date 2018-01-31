#!/usr/bin/python
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from analyze_linear import trapezoid, measure_error
from analyze_amr_linear import read_amr

L2, x, fixed, analytic = measure_error("linear.640.txt")

x_amr, amr = read_amr(["linear_amr.10.0.txt","linear_amr.20.0.txt","linear_amr.40.0.txt","linear_amr.40.7.txt","linear_amr.40.14.txt","linear_amr.80.14.txt","linear_amr.80.28.txt","linear_amr.160.27.txt","linear_amr.160.54.txt","linear_amr.320.54.txt","linear_amr.320.108.txt","linear_amr.640.107.txt","linear_amr.640.214.txt"])

x_amr = np.array(x_amr)

amr_analytic = np.zeros(len(amr))
amr_analytic[np.where(np.array(x_amr)<0.75)] = 1.0

amr_L2 = trapezoid(x_amr, (amr - amr_analytic)**2)
 
print "fixed", L2
print "amr", amr_L2

plt.figure()
plt.ylabel("phi")
plt.xlabel("x")
plt.ylim([-0.1,1.1])
plt.xlim([-0.1,1.1])
plt.title("Level 7 Comparison")
plt.plot(x, fixed, 'x', markersize=15, label="Fixed", color='g')
plt.plot(x_amr, amr, '.', markersize=15, label="AMR", color='r')
plt.plot(x, analytic, linewidth=2, label="analytic")
plt.legend(loc='best')

font = {'family' : 'normal',
        'weight' : 'bold',
        'size' : 18}

matplotlib.rc('font',**font)
plt.show()

