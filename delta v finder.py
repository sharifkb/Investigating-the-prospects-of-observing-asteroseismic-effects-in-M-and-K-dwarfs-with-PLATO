# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 04:18:11 2020

@author: shari
"""

import numpy as np
import  astropy.constants as c

#constants
L_sun = c.L_sun.value
R_sun = c.R_sun.to("km").value 
delta_v_sun = 134.8
sigma1 = c.sigma_sb.value

"""
Change this to your relevant directory
"""
###### M DWARFS
#catalogue = np.load(r"C:\Users\shari\OneDrive\Documents\Work\3rd Yr\Group studies\Code\Data array creating code\Mdwarf\M_dwarves_plat_1_final.npy")

###### K DWARFS
catalogue = np.load(r"C:\Users\shari\OneDrive\Documents\Work\3rd Yr\Group studies\Code\Data array creating code\Kdwarf\K_dwarves_plat_1_final.npy")

logTe = catalogue["logTe"]
logL = catalogue["logL"]
M = catalogue["Mass"]
def spectrum (logTe, M, logL):
    
    T_eff = 10**logTe
    L1 = 10**logL*L_sun


    R = np.sqrt(L1/(4*np.pi*sigma1*T_eff**4))*10**(-3)
    delta_v = delta_v_sun*M**(1/2)*(R/R_sun)**(-3/2)

    return(delta_v)

"""
Creating a catalogue of multiple stars. Change the number of stars you want to catalogue
"""
numstars = 1128
start = 0
num = range(1128)
fish = num[start:start+numstars]

V_max = np.empty((len(fish)))


for i in fish:

    V_max[i-start] = spectrum(logTe[i], M[i], logL[i])


dat1 = np.vstack([V_max.T]).T
###219.029884120290
###145.77899479508767
"""
Change to your relevant directory
"""
###### M DWARFS
#np.savetxt(r"C:\Users\shari\OneDrive\Documents\Work\3rd Yr\Group studies\Code\Data array creating code\Mdwarf\Mypoints"+str(fish)+".txt", dat0, delimiter = ' ')
#np.savetxt(r"C:\Users\shari\OneDrive\Documents\Work\3rd Yr\Group studies\Code\Data array creating code\Mdwarf\Motherstuff"+str(fish)+".txt", dat1, delimiter = ' ', fmt='%s' )
###### K DWARFS
#np.savetxt(r"C:\Users\shari\OneDrive\Documents\Work\3rd Yr\Group studies\Code\Data array creating code\Kdwarf\Kypoints"+str(fish)+".txt", dat0, delimiter = ' ')
#np.savetxt(r"C:\Users\shari\OneDrive\Documents\Work\3rd Yr\Group studies\Code\Data array creating code\Kdwarf\Kotherstuff"+str(fish)+".txt", dat1, delimiter = ' ', fmt='%s' )
