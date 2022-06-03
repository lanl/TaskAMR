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
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from analyze_linear import trapezoid, measure_error
from analyze_amr_linear import read_amr

L2, x, fixed, analytic = measure_error("linear.80.txt")

filenames = ["0319_linear_amr.20.0.txt","0319_linear_amr.40.0.txt", "0319_linear_amr.80.6.txt","0319_linear_amr.40.3.txt","0319_linear_amr.80.12.txt", "0319_linear_amr.80.18.txt","0319_linear_amr.80.24.txt"]
x_amr, amr = read_amr(filenames)
print(x_amr)
x_amr = np.array(x_amr)

amr_analytic = np.zeros(len(amr))
amr_analytic[np.where(np.array(x_amr)<0.75)] = 1.0

amr_L2 = trapezoid(x_amr, (amr - amr_analytic)**2)
 
print("fixed", L2)
print("amr", amr_L2)

plt.figure()
plt.ylabel("phi")
plt.xlabel("x")
plt.ylim([-0.1,1.1])
plt.xlim([-0.1,1.1])
plt.title("Level 4 Comparison")
plt.plot(x, fixed, 'x', markersize=15, label="Fixed", color='g')
plt.plot(x_amr, amr, '.', markersize=15, label="AMR", color='r')
plt.plot(x, analytic, linewidth=2, label="analytic")
plt.legend(loc='best')

font = {'family' : 'normal',
        'weight' : 'bold',
        'size' : 18}

matplotlib.rc('font',**font)
plt.show()

