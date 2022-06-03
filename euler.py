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

GAMMA = 1.4
BETA = (GAMMA + 1.0) / (GAMMA - 1.0)


x = np.arange(0.0, 1.01, .01)
x_jump = 0.5
t_final = 0.142681382

# left and right states
P_l = 1.0
rho_l = 1.0
u_l = 0.0

P_r = 0.1
rho_r = 0.125
u_r = 0.0

P_0 = np.zeros(len(x))
rho_0 = np.zeros(len(x))
u_0 = np.zeros(len(x))

x_l = np.where(x < 0.5)
x_r = np.where(x >= 0.5)

P_0[x_l] = P_l
rho_0[x_l] = rho_l
u_0[x_l] = u_l

P_0[x_r] = P_r
rho_0[x_r] = rho_r
u_0[x_r] = u_r

def plot_density(density):
    import matplotlib.pyplot as plt
    plt.figure()
    plt.ylabel("density")
    plt.yticks([0,0.3,0.6,0.9,1.2])
    plt.ylim([0,1.2])
    plt.xlabel("x")
    plt.xticks([0,0.25,0.5,0.75,1.])
    plt.plot(x,density)

def plot_speed(speed):
    import matplotlib.pyplot as plt
    plt.figure()
    plt.ylabel("speed")
    #plt.yticks([0,0.3,0.6,0.9,1.2])
    #plt.ylim([0,1.2])
    plt.xlabel("x")
    plt.xticks([0,0.25,0.5,0.75,1.])
    plt.plot(x,speed)


def plot_velocity(velocity):
    import matplotlib.pyplot as plt
    plt.figure()
    plt.ylabel("velocity")
    plt.yticks([0,0.3,0.6,0.9,1.2])
    plt.ylim([0,1.2])
    plt.xlabel("x")
    plt.xticks([0,0.25,0.5,0.75,1.])
    plt.plot(x,velocity)


def plot_pressure(pressure):
    import matplotlib.pyplot as plt
    plt.figure()
    plt.ylabel("pressure")
    plt.yticks([0,0.3,0.6,0.9,1.2])
    plt.ylim([0,1.2])
    plt.xlabel("x")
    plt.xticks([0,0.25,0.5,0.75,1.])
    plt.plot(x,pressure)


def plot_sie(sie):
    import matplotlib.pyplot as plt
    plt.figure()
    plt.ylabel("specific internal energy")
    plt.yticks([1,1.5,2.0,2.5,3.0])
    plt.ylim([1,3.0])
    plt.xlabel("x")
    plt.xticks([0,0.25,0.5,0.75,1.])
    plt.plot(x,sie)

def specific_internal_energy(P,rho):
    return P / (rho * (GAMMA - 1.0))

def get_energy(P, rho, v):
    return P / (GAMMA - 1.0) + 0.5 * rho * v**2

def get_pressure(E, rho, v):
    return (E - 0.5 * rho * v**2) * (GAMMA-1.0)

def get_flux(E, rho, momentum):
    v = momentum / rho
    value = np.zeros((3, len(E)))
    value[0,:] = rho * v
    P = get_pressure(E, rho, v)
    value[1,:] = rho * v**2 + P
    value[2,:] = (E + P) * v
    return value

