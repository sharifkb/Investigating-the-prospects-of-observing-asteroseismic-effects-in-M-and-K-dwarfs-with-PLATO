# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 21:05:12 2020

@author: shari
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import  astropy.constants as c

data = pd.read_table("synthetic_trilegal.dat", sep="\s+")
df = pd.DataFrame(data)

B = df.iloc[:, 12]
V = df.iloc[:, 13]
L = df.iloc[:, 4]
mbol = df.iloc[:, 10]
distmod = df.iloc[:, 7]
T_eff = df.iloc[:, 5]

M_bol11 = -2.5*(L+np.log10(c.L_sun/c.L_bol0))
M_bol1 = -2.5*L
M_bol2 = mbol-distmod
M_bol_err = np.sqrt((M_bol1-M_bol2)**2)
Col = B-V

plt.figure()
plt.plot(Col, M_bol1, "r.")
#plt.plot(Col, M_bol11, "y.")
plt.gca().invert_yaxis()
plt.grid()
plt.title("Colour Magnitude Diagram")
plt.ylabel("Absolute Magnitude")
plt.xlabel("B-V Colour")

plt.figure()
plt.plot(Col, mbol, "r.")
plt.gca().invert_yaxis()
plt.grid()
plt.title("Colour Magnitude Diagram")
plt.ylabel("Apparent Bolometric Magnitude")
plt.xlabel("B-V Colour")

plt.figure()
axes = plt.gca()
plt.plot(T_eff, L, "r.")
plt.gca().invert_xaxis()
plt.grid()
axes.set_xlim([4.6,3.4])
plt.title("Hertzsprung Russell Diagram")
plt.ylabel("Luminosity $log(L/L_{\odot})$")
plt.xlabel("Effective Temperature $log(T_{eff})$")
