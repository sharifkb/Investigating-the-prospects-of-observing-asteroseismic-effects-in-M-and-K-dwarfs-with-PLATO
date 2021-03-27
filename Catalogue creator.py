# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 17:11:24 2020

@author: shari
"""
import numpy as np
import  astropy.constants as c

#constants
tau_sun = 200e-6
g_sun = 274
T_eff_sun = 5777
sigma_sun = 50000/(c.R_sun.value/(1000*10**3))
L_sun = c.L_sun.value
M_sun = c.M_sun.value #kg
R_sun = c.R_sun.to("km").value 
V_max_sun = 3100
delta_v_sun = 134.8
age_sun = 4.569*10**9
A_sun = 2.1
D_sun = 1.493
n = np.linspace(1,40,40)
l0= 0
l1 = 1
l2 = 2
l3 = 3
epsilon=1.4426
a_0 = -25.5
b_0 = 29.1
a_2 = 6.3
b_2 = -1.8
sigma1 = c.sigma_sb.value
size = 20000
space = size*2*365*24*3600/1e6
x2 = np.linspace(1, size, space)
sum_vis = 0.5 + 1.5 + 0.04

"""
Change this to your relevant directory
"""
###### M DWARFS
catalogue = np.load(r"C:\Users\shari\OneDrive\Documents\Work\3rd Yr\Group studies\Code\Data array creating code\Mdwarf\M_dwarves_plat_1_final.npy")

###### K DWARFS
#catalogue = np.load(r"C:\Users\shari\OneDrive\Documents\Work\3rd Yr\Group studies\Code\Data array creating code\Kdwarf\K_dwarves_plat_1_final.npy")

logg = catalogue["logg"]
logTe = catalogue["logTe"]
logL = catalogue["logL"]
M = catalogue["Mass"]
vmag =  catalogue["V"]
age = catalogue["logAge"]

def spectrum (logTe, logg, M, logL, vmag):
    
    T_eff = 10**logTe
    g = 10**logg/100
    L1 = 10**logL*L_sun
    L = 10**logL

    R = np.sqrt(L1/(4*np.pi*sigma1*T_eff**4))*10**(-3)
    delta_v = delta_v_sun*M**(1/2)*(R/R_sun)**(-3/2)
    D = D_sun*(delta_v/delta_v_sun)
    V_max = V_max_sun*(g/g_sun)*(T_eff/T_eff_sun)**(-1/2)
    T_red = 8907*L**-0.093
    beta = 1 - np.exp(-(T_red - T_eff)/1250)
    A_max = beta*A_sun*L/M*(T_eff/T_eff_sun)**(-2)
    sigma = 2.123*delta_v
    Gauss = A_max**2*np.exp(-0.5*((x2-V_max)/sigma)**2)
    
    linewidth = 0.07+0.91*(T_eff/5777)**15.3
    height = (Gauss*2*delta_v)/(linewidth*np.pi*sum_vis)

    tau_c = tau_sun*(g_sun/g)*(T_eff/T_eff_sun)**(1/2)
    sigma_c = sigma_sun*(T_eff/T_eff_sun)*R/(M*R_sun)
    gran = 4*np.sqrt(2)*tau_c*sigma_c**2/(1+(2*np.pi*x2*tau_c)**4)             
    
    vmag_value = np.full(len(x2),vmag)
    shotnoise = 18.0*10**(-0.4*(11.0-vmag_value))
    
    modes0=[]
    modes1=[]
    modes2=[]
    modes3=[]
    for i in range(len(n)):
        modes0.append(delta_v*(n[i]+l0/2+epsilon)-D*l0*(l0+1))
        modes1.append(delta_v*(n[i]+l1/2+epsilon)-D*l1*(l1+1))  
        modes2.append(delta_v*(n[i]+l2/2+epsilon)-D*l2*(l2+1))  
        modes3.append(delta_v*(n[i]+l3/2+epsilon)-D*l3*(l3+1))    
    
    a0 = np.zeros(len(x2))
    a1 = np.zeros(len(x2))
    a2 = np.zeros(len(x2))
    a3 = np.zeros(len(x2))
    for i in range(0,40):
        a0+=height/(1+((x2 - np.array(modes0[i]))/(linewidth/2))**2)
        a1+=height*1.5/(1+((x2 - np.array(modes1[i]))/(linewidth/2))**2)   
        a2+=height*0.5/(1+((x2 - np.array(modes2[i]))/(linewidth/2))**2)   
        a3+=height*0.04/(1+((x2 - np.array(modes3[i]))/(linewidth/2))**2)               

    signal0 = a0 + +a1 + a2 +a3 + gran + shotnoise

    return(signal0)

def otherstuff (logTe, logg, M, logL, vmag):
    
    T_eff = 10**logTe
    g = 10**logg/100
    L1 = 10**logL*L_sun

    R = np.sqrt(L1/(4*np.pi*sigma1*T_eff**4))*10**(-3)
    V_max = V_max_sun*(g/g_sun)*(T_eff/T_eff_sun)**(-1/2)
    delta_v = delta_v_sun*M**(1/2)*(R/R_sun)**(-3/2) 
    sigma = 2.123*delta_v
    return(V_max, sigma)
    

"""
Creating a catalogue of multiple stars. Change the number of stars you want to catalogue
"""
numstars = 1
start = 0
num = range(1128)
fish = num[start:start+numstars]

star = np.empty((len(x2), len(fish)))
V_max = np.empty((len(fish)))
vmagn = np.empty((len(fish)))
width = np.empty((len(fish)))
Temp = np.empty((len(fish)))

for i in fish:
    star[:,i-start] = spectrum(logTe[i], logg[i], M[i], logL[i], vmag[i])
    V_max[i-start] = otherstuff(logTe[i], logg[i], M[i], logL[i], vmag[i])[0]
    width[i-start] = otherstuff(logTe[i], logg[i], M[i], logL[i], vmag[i])[1]
    vmagn[i-start] = vmag[i]
    Temp[i-start] = logTe[i]
dat0 = np.vstack([x2, star.T])
dat1 = np.vstack([fish, V_max.T, vmagn, width, Temp]).T

"""
Change to your relevant directory
"""
###### M DWARFS
#np.savetxt(r"C:\Users\shari\OneDrive\Documents\Work\3rd Yr\Group studies\Code\Data array creating code\Mdwarf\Mypoints"+str(fish)+".txt", dat0, delimiter = ' ')
#np.savetxt(r"C:\Users\shari\OneDrive\Documents\Work\3rd Yr\Group studies\Code\Data array creating code\Mdwarf\Motherstuff"+str(fish)+".txt", dat1, delimiter = ' ', fmt='%s' )
###### K DWARFS
#np.savetxt(r"C:\Users\shari\OneDrive\Documents\Work\3rd Yr\Group studies\Code\Data array creating code\Kdwarf\Kypoints"+str(fish)+".txt", dat0, delimiter = ' ')
#np.savetxt(r"C:\Users\shari\OneDrive\Documents\Work\3rd Yr\Group studies\Code\Data array creating code\Kdwarf\Kotherstuff"+str(fish)+".txt", dat1, delimiter = ' ', fmt='%s' )
"""
For cataloging one star at a time. Change the value of i to choose the star you want.
"""
