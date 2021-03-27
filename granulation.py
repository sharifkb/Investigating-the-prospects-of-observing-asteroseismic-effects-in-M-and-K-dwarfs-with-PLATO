# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 18:43:25 2020

@author: shari
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import  astropy.constants as c
import scipy.constants as d
from tomso import fgong
import matplotlib as mpl
mpl.rc('font', family='serif')

#https://users-phys.au.dk/~jcd/solar_models/
#https://tomso.readthedocs.io/en/stable/fgong.html

S = fgong.load_fgong(r"C:\Users\shari\OneDrive\Documents\Work\3rd Yr\Group studies\Code\fgong.l5bi.d.15", return_object=True)
catalogue = np.load(r"C:\Users\shari\OneDrive\Documents\Work\3rd Yr\Group studies\Code\Data array creating code\Kdwarf\K_dwarves_plat_1_final.npy")

#stellar values in cgs
r      = S.r       # Radius
r2     = S.x
s      = S.cs      # speed of sound
rho    = S.rho     # density
p      = S.P       # pressure
T      = S.T       # Temperature
k      = S.kappa   # Opacity
L      = S.L_r     # Luminosity as a function of r
m      = S.m       # Mass as a function of r   
gamma  = S.G1      # First adiabatic index as a function of r
H      = S.X       # Fractional hydrogen abundance as a function of r

#### Assuming hydorgen and helium only
X  = S.X
Z  = S.Z
Y  = 1-X-Z
mu = (2*X+3*Y/4+Z/2)**-1
#(H+He*He_m)

# Constants have been converted to cgs by *10**n

nabla_ad  = (gamma-1)/gamma
nabla     = np.gradient(np.log(T), np.log(p))
nabla_mu  = np.gradient(np.log(mu), np.log(p))
a         = 4*d.sigma*10**3    #a factor of c had been omitted since it would've been multiplied out in the next line
C         = 3*d.k*10**7/(16*d.pi*a*c.G.cgs.value*d.m_p*10**3)
nabla_rad = C*k*L*rho/(mu*m*T**3)

logTe = catalogue["logTe"]
logL  = catalogue["logL"]
Mass  = catalogue["Mass"]

n = 20
def spectrum (logTe, M, logL):
     
    T_eff = 10**logTe
    M     = M*c.M_sun.cgs.value
    L    = 10**logL*c.L_sun.cgs.value  #values are given in multiples of solar luminosities
    R     = (L/(4*np.pi*d.sigma*10**3*T_eff**4))**(1/2)
    #scale = (M/c.M_sun.cgs.value)*(R/c.R_sun.cgs.value)**-4*(T_eff/5777)**(-5/2)  
    #scale = 1
    scale = (5777/T_eff)**(7/2)
    print(T_eff)
    
    # All temperature gradients are scaled from the sun using this scaling
    y1 = nabla
    y2 = nabla_ad
    y3 = nabla_rad*scale
    y4 = nabla_mu
    
    fig1  = plt.figure(figsize=(8,5.5))
    #plt.plot(r2, y1, label = r"$\nabla$", color = "green")
    plt.plot(r2, y2, label = r"$\nabla_{ad}$", color = "red")
    plt.plot(r2, y3, label = r"$\nabla_{rad}$", color = "blue")
    idx = np.argwhere(np.diff(np.sign(y2 - y3))).flatten()
    x1 = r2[idx][0]
    x2 = r2[idx][1]
    plt.axvline(x = x1, color = "k", linestyle = "--")
    plt.axvline(x = x2, color = "k", linestyle = "--")
    plt.axvspan(x2, x1, alpha=0.5, label = 'Convective zone', color='grey')
    
    plt.ylim(0, np.max(nabla_ad*scale)*2)
    plt.ylabel(r"$\nabla$", fontsize = n)
    plt.xlabel(r"r/R", fontsize = n)
    plt.tick_params(axis='both', which='major', labelsize=n)
    plt.legend(fontsize = n) 
    plt.tight_layout()

    return(fig1)

numstars = 1
start    = 6
num      = range(1128)
fish     = num[start:start+numstars]

for i in fish:
    spectrum(logTe[i], Mass[i], logL[i])
#sun
#spectrum(3.76133, 1, 0)

