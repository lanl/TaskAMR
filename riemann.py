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

from euler import GAMMA, BETA, x, x_jump, t_final, P_l, rho_l, u_l, P_r, rho_r, u_r
from euler import plot_density, plot_velocity, plot_pressure, plot_sie, plot_speed
from euler import specific_internal_energy, get_energy

EPS = 1.0e-6

def speed_of_sound(P, rho):
    return np.sqrt(GAMMA * P / rho)

def phi_rarefact(P_star, P, rho):
    c = speed_of_sound(P, rho)
    phi = 2.0*c*(1.0-(P_star/P)**(0.5*(GAMMA-1.0)/GAMMA))/(GAMMA-1.0) # (14.51) Leveque
    phi[np.where(P_star>P)] = np.nan
    return phi

def phi_shock(P_star, P, rho):
    c = speed_of_sound(P, rho)
    phi = 2.0*c * (1.0 - P_star/P) / np.sqrt(2.0*GAMMA*(GAMMA-1)*(1.0+BETA*P_star/P)) # (14.52) Leveque
    phi[np.where(P_star<P)] = np.nan
    return phi

def one_rarefact(P, P_l, rho_l, u_l):
    u = u_l + phi_rarefact(P, P_l, rho_l) # (14.51) Leveque
    return u

def one_shock(P, P_l, rho_l, u_l):
    u = u_l + phi_shock(P, P_l, rho_l) # (14.52) Leveque
    return u

def three_rarefact(P, P_r, rho_r, u_r):
    u = u_r - phi_rarefact(P, P_r, rho_r) # (14.54) Leveque
    return u

def three_shock(P, P_r, rho_r, u_r):
    u = u_r - phi_shock(P, P_r, rho_r) # (14.55) Leveque
    return u

def phi(P_star, P, rho):
    val = np.nan
    if (P_star>=P):
        val = phi_shock(np.array([P_star]), P, rho)
    else:
        val = phi_rarefact(np.array([P_star]), P, rho)
    return val[0]

def one(P, P_l, rho_l, u_l):
    return u_l + phi(P, P_l, rho_l)

def three(P, P_r, rho_r, u_r):
    return u_r - phi(P, P_r, rho_r)

def deriv_phi_rarefact(P_star, P, rho):
    c = speed_of_sound(P, rho)
    val = -c * (P_star/P)**(-0.5*(GAMMA+1.0)/GAMMA) /GAMMA
    return val

def deriv_phi_shock(P_star, P, rho):
    c = speed_of_sound(P, rho)
    val = -2.0*c * (1.0 + 0.5 * BETA * (1.0-P_star/P) / (1.0+BETA*P_star/P) ) / (np.sqrt(2.0*GAMMA*(GAMMA-1)*(1.0+BETA*P_star/P)) * P)
    return val

def deriv_phi(P_star, P, rho):
    val = np.nan
    if (P_star>=P):
        val = deriv_phi_shock(P_star, P, rho)
    else:
        val = deriv_phi_rarefact(P_star, P, rho)
    return val

def rho_star(P_star, P, rho):
    val = np.nan
    if (P_star>=P):
        val = (1.0+BETA*P_star/P)*rho/(P_star/P + BETA)
    else:
        val = (P_star / P)**(1.0/GAMMA) * rho
    return val

def verify_Rankine_Hugoniot(P_u, rho_u, v_u, P_d, rho_d, v_d, S):
    E_u = get_energy(P_u, rho_u, v_u)
    E_d = get_energy(P_d, rho_d, v_d)
    #print(v_u * (E_u + P_u) - v_d*(E_d+P_d),"vs",S*(E_u - E_d))
    #print(rho_u*v_u**2+P_u - rho_d*v_d**2 - P_d,"vs",S*(rho_u*v_u-rho_d*v_d))
    #print(rho_u*v_u - rho_d*v_d,"vs",S*(rho_u-rho_d))

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    P = np.arange(0.01, 2.0, 0.01)

    plt.figure()
    plt.ylabel("u = velocity")
    plt.xlabel("P = pressure")
    plt.xlim([0,2])
    plt.plot([P_l],[u_l],'+',label='(P_l,u_l)')
    plt.plot([P_r],[u_r],'x',label='(P_r,u_r)')
    plt.plot(P,one_shock(P, P_l, rho_l, u_l),"--",label='1-shock')
    plt.plot(P,one_rarefact(P, P_l, rho_l, u_l),":",label='1-rarefact')
    plt.plot(P,three_shock(P, P_r, rho_r, u_r),"--",label='3-shock')
    plt.plot(P,three_rarefact(P, P_r, rho_r, u_r),":",label='3-rarefact')
    plt.legend(loc='best')

    #plt.figure()
    #plt.ylabel("f = u_l^* - u_r^*")
    #plt.xlabel("P = pressure")
    #plt.xlim([0,2])
    #f_func = np.zeros(len(P))
    #for i in range(len(P)):
        #f_func[i] = one(P[i], P_l, rho_l, u_l) - three(P[i], P_r, rho_r, u_r)
    #plt.plot(P,f_func)

    # Newton solve
    P_star = 0.5 * (P_l + P_r)
    u_l_star = one(P_star, P_l, rho_l, u_l)
    u_r_star = three(P_star, P_r, rho_r, u_r)
    f = u_l_star - u_r_star

    count = 0
    while(np.abs(f) > EPS):
        #print('du_l',deriv_phi(P_star, P_l, rho_l))
        #print('du_r', -1.0 * deriv_phi(P_star, P_r, rho_r))
        # u_l_star = u_l + phi, u_r_star = u_r - phi
        df = deriv_phi(P_star, P_l, rho_l) + deriv_phi(P_star, P_r, rho_r)
        print('iter',count,'f',f,'df',df)

        # Newton update
        P_star = P_star - f/df
        count += 1

        u_l_star = one(P_star, P_l, rho_l, u_l)
        u_r_star = three(P_star, P_r, rho_r, u_r)
        plt.plot([P_star],[u_l_star], 'o')
        plt.plot([P_star],[u_r_star], '.')
        f = u_l_star - u_r_star

    print('P_star',P_star,'u_l_star',u_l_star,'u_r_star',u_r_star)
    u_star = 0.5 * (u_l_star + u_r_star)

    S_3 = (P_star - P_r - rho_r * u_r**2)/(rho_r * (u_star - u_r))
    print("S_3", S_3, "u+c",u_r + speed_of_sound(P_r, rho_r))


    rho_l_star = rho_star(P_star, P_l, rho_l)
    print('rho_l_star',rho_l_star)
    rho_r_star = rho_star(P_star, P_r, rho_r)
    print('rho_r_star',rho_r_star,'or', S_3 * rho_r / (S_3 - u_star))

    print
    verify_Rankine_Hugoniot(P_r, rho_r, u_r, P_star, rho_r_star, u_star, S_3)
    print

    # Sod's time
    delta_t = t_final

    S_2 = u_star
    print("S_2", S_2)
    c_l = speed_of_sound(P_l, rho_l)
    S_1_head = u_l - c_l
    print("S_1_head",S_1_head)
    S_1_tail = u_star - speed_of_sound(P_star, rho_l_star)
    print("S_1_tail",S_1_tail)

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

    plot_density(density)
    plot_velocity(velocity)
    plot_pressure(pressure)
    plot_sie(sie)
    plot_speed(speed_of_sound(pressure,density))

    plt.show()
