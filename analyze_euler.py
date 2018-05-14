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

from euler import GAMMA, BETA, x, x_jump, P_l, rho_l, u_l, P_r, rho_r, u_r
from riemann import one, three, EPS, deriv_phi, rho_star, verify_Rankine_Hugoniot
from riemann import speed_of_sound

t_final = 0.142625

def specific_internal_energy(P, rho):
  return P / (rho * (GAMMA - 1.0))

def get_pressure(E, rho, v):
  return (E - 0.5 * rho * v**2) * (GAMMA-1.0)

def reimann_solve(x):
    # Newton solve
    P_star = 0.5 * (P_l + P_r)
    u_l_star = one(P_star, P_l, rho_l, u_l)
    u_r_star = three(P_star, P_r, rho_r, u_r)
    f = u_l_star - u_r_star
    count = 0
    while(np.abs(f) > EPS):
      df = deriv_phi(P_star, P_l, rho_l) + deriv_phi(P_star, P_r, rho_r)
      # Newton update
      P_star = P_star - f/df
      count += 1
      u_l_star = one(P_star, P_l, rho_l, u_l)
      u_r_star = three(P_star, P_r, rho_r, u_r)
      f = u_l_star - u_r_star
    u_star = 0.5 * (u_l_star + u_r_star)
    S_3 = (P_star - P_r - rho_r * u_r**2)/(rho_r * (u_star - u_r))
    rho_l_star = rho_star(P_star, P_l, rho_l)
    rho_r_star = rho_star(P_star, P_r, rho_r)
    verify_Rankine_Hugoniot(P_r, rho_r, u_r, P_star, rho_r_star, u_star, S_3)
    # Sod's time
    delta_t = t_final
    S_2 = u_star
    c_l = speed_of_sound(P_l, rho_l)
    S_1_head = u_l - c_l
    S_1_tail = u_star - speed_of_sound(P_star, rho_l_star)
    density = np.zeros(len(x))
    velocity = np.zeros(len(x))
    pressure = np.zeros(len(x))
    sie = np.zeros(len(x))
    for i in range(len(x)):
      if ((x[i]-x_jump) <= (delta_t*S_1_head)):
        density[i] = rho_l
        velocity[i] = u_l
        pressure[i] = P_l
      elif ((x[i]-x_jump) > (delta_t*S_1_head)) & ((x[i]-x_jump)<=(delta_t*S_1_tail)):
        Xsi = (x[i]-x_jump) / delta_t
        velocity[i] = ((GAMMA-1.0)*u_l+2.0*(c_l+Xsi)) / (GAMMA+1.0)
        density[i] = (rho_l**GAMMA*(velocity[i]-Xsi)**2/(GAMMA*rho_l))**(1.0/(GAMMA-1.0))
        pressure[i] = P_l/rho_l**GAMMA * density[i]**GAMMA
      elif ((x[i]-x_jump) > (delta_t*S_1_tail)) & ((x[i]-x_jump)<=(delta_t*S_2)):
        density[i] = rho_l_star
        velocity[i] = u_star
        pressure[i] = P_star
      elif ((x[i]-x_jump) > (delta_t*S_2)) & ((x[i]-x_jump)<=(delta_t*S_3)):
        density[i] = rho_r_star
        velocity[i] = u_star
        pressure[i] = P_star
      elif (x[i]-x_jump) > (delta_t*S_3):
        density[i] = rho_r
        velocity[i] = u_r
        pressure[i] = P_r
      sie[i] = specific_internal_energy(pressure[i], density[i])
    return density, velocity, pressure, sie

def measure_error(filename):
    with open(filename,"r") as f:
      density = []
      momentum = []
      energy = []
      for line in f:
        data = line.split()
        density.append(float(data[0]))
        momentum.append(float(data[1]))
        energy.append(float(data[2]))
      x = (0.5 + np.arange(float(len(density))) )/float(len(density))

      momentum = np.array(momentum)
      num_density = np.array(density)
      energy = np.array(energy)
      num_velocity = momentum / num_density
      num_pressure = get_pressure(energy, num_density, num_velocity)
      num_sie = specific_internal_energy(num_pressure, num_density)

      density, velocity, pressure, sie = reimann_solve(x)

      sie_L2 = np.mean((num_sie - sie)**2)
      P_L2 = np.mean((num_pressure - pressure)**2)
      v_L2 = np.mean((num_velocity - velocity)**2)
      rho_L2 = np.mean((num_density - density)**2)
      L2 = np.mean([sie_L2, P_L2, v_L2, rho_L2])

    return L2, x, num_density, density

if __name__== "__main__":
 
  parser = argparse.ArgumentParser(description='Plot convergence for fixed grid linear advection.')
  parser.add_argument('text_files',nargs='*')

  args = parser.parse_args()
  NX = []
  Error = []

  plt.figure()
  plt.ylabel("density")
  plt.yticks([0,0.3,0.6,0.9,1.2])
  plt.ylim([0,1.2])
  plt.xlabel("x")
  plt.xticks([0,0.25,0.5,0.75,1.])
  x = 1.0e-3 * (0.5 + np.arange(float(1000)) )
  density, velocity, pressure, sie = reimann_solve(x)
  plt.plot(x,density,label='anal')

  for filename in args.text_files:
    print filename
    resolution = filename.split('.')
    NX.append(float(resolution[1]))
    L2, x, num_density, density = measure_error(filename)
    plt.plot(x,num_density,'--',label='NX='+str(NX[-1]))
    Error.append(L2)

  print NX
  print Error
  plt.legend(loc='best')

  plt.figure()
  plt.title("Convergence plot")
  plt.ylabel("L_2 error")
  plt.xlabel("NX")
  plt.yscale('log')
  plt.xscale('log')
  plt.plot(NX,Error,'o')

  fit = np.polyfit(np.log(NX[5:]),np.log(Error[5:]),1)
  p = fit[0]
  A = np.exp(fit[1])
  plt.plot(NX, A * NX**p,label="NX^"+str(p))
  plt.legend(loc='best')
  plt.show()

