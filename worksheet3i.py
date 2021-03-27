# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 15:04:26 2020

@author: shari
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

num = [1,2,3,4,5]

#names = ["l", "n", "freq", "A"]
s1 = pd.read_table("m1.00Y0.278Z0.02a1.90.log3.freq", sep="\s+", skiprows=2)
s2 = pd.read_table("m1.00Y0.278Z0.02a1.90.log5.freq", sep="\s+", skiprows=2)
s3 = pd.read_table("m1.00Y0.278Z0.02a1.90.log7.freq", sep="\s+", skiprows=2)
s4 = pd.read_table("m1.00Y0.278Z0.02a1.90.log9.freq", sep="\s+", skiprows=2)
s5 = pd.read_table("m1.00Y0.278Z0.02a1.90.log11.freq", sep="\s+", skiprows=2)

df1 = pd.DataFrame(s1)
df2 = pd.DataFrame(s2)
df3 = pd.DataFrame(s3)
df4 = pd.DataFrame(s4)
df5 = pd.DataFrame(s5)

a = df1.iloc[:, 2]
b = df2.iloc[:, 2]
c = df3.iloc[:, 2]
d = df4.iloc[:, 2]
e = df5.iloc[:, 2]

mass_frac = [np.log10(6.93e-01), np.log10(4.42e-01), np.log10(2.23e-01), np.log10(1.05e-03), np.log10(7.25e-07)]


def av (a, b):
    return(a+b)/2

def diff(a, b):
    return np.sqrt((a-b)**2)

##########################################################################################
#star 1
afreq1 = np.array(a)[0:13]
afreq2 = np.array(a)[1:14]
afreq3 = np.array(a)[15:28]
afreq4 = np.array(a)[16:29]
afreq5 = np.array(a)[30:43]
afreq6 = np.array(a)[31:44]


asep1 = diff(afreq2, afreq1)
asep2 = diff(afreq4, afreq3)
asep3 = diff(afreq6, afreq5)
asep4 = av(afreq2, afreq1)
asep5 = av(afreq4, afreq3)
asep6 = av(afreq6, afreq5)

adeltfreq1 = np.array(a)[0:14]
adeltfreq5 = np.array(a)[30:44]

adeltsep1 = diff(adeltfreq1, adeltfreq5)
adeltsep2 = av(adeltfreq1, adeltfreq5)


plt.figure()
plt.plot(asep4, asep1, "r", label="Angular degree 0")
plt.plot(asep5, asep2, "g", label="Angular degree 1")
plt.plot(asep6, asep3, "b", label="Angular degree 2")
plt.title("Log3 Large Separation as a function of Frequency")
plt.ylabel("Large Separation, $\Delta$v, $\mu$Hz")
plt.xlabel("Frequency, $\mu$Hz")
plt.legend(loc="upper left")

##########################################################################################
#star 2
bfreq1 = np.array(b)[0:13]
bfreq2 = np.array(b)[1:14]
bfreq3 = np.array(b)[15:27]
bfreq4 = np.array(b)[16:28]
bfreq5 = np.array(b)[29:42]
bfreq6 = np.array(b)[30:43]

bsep1 = diff(bfreq2, bfreq1)
bsep2 = diff(bfreq4, bfreq3)
bsep3 = diff(bfreq6, bfreq5)
bsep4 = av(bfreq2, bfreq1)
bsep5 = av(bfreq4, bfreq3)
bsep6 = av(bfreq6, bfreq5)

bdeltfreq1 = np.array(b)[0:14]
bdeltfreq5 = np.array(b)[29:43]

bdeltsep1 = diff(bdeltfreq1, bdeltfreq5)
bdeltsep2 = av(bdeltfreq1, bdeltfreq5)

plt.figure()
plt.plot(bsep4, bsep1, "r", label="Angular degree 0")
plt.plot(bsep5, bsep2, "g", label="Angular degree 1")
plt.plot(bsep6, bsep3, "b", label="Angular degree 2")
plt.title("Log5 Large Separation as a function of Frequency")
plt.ylabel("Large Separation, $\Delta$v, $\mu$Hz")
plt.xlabel("Frequency, $\mu$Hz")
plt.legend(loc="upper left")

##########################################################################################
#star 3
cfreq1 = np.array(c)[0:12]
cfreq2 = np.array(c)[1:13]
cfreq3 = np.array(c)[14:26]
cfreq4 = np.array(c)[15:27]
cfreq5 = np.array(c)[28:39]
cfreq6 = np.array(c)[29:40]

csep1 = diff(cfreq2, cfreq1)
csep2 = diff(cfreq4, cfreq3)
csep3 = diff(cfreq6, cfreq5)
csep4 = av(cfreq2, cfreq1)
csep5 = av(cfreq4, cfreq3)
csep6 = av(cfreq6, cfreq5)

cdeltfreq1 = np.array(c)[1:13]
cdeltfreq5 = np.array(c)[28:40]

cdeltsep1 = diff(cdeltfreq1, cdeltfreq5)
cdeltsep2 = av(cdeltfreq1, cdeltfreq5)

plt.figure()
plt.plot(csep4, csep1, "r", label="Angular degree 0")
plt.plot(csep5, csep2, "g", label="Angular degree 1")
plt.plot(csep6, csep3, "b", label="Angular degree 2")
plt.title("Log7 Large Separation as a function of Frequency")
plt.ylabel("Large Separation, $\Delta$v, $\mu$Hz")
plt.xlabel("Frequency, $\mu$Hz")
plt.legend(loc="upper left")

##########################################################################################
#star 4
dfreq1 = np.array(d)[0:11]
dfreq2 = np.array(d)[1:12]
dfreq3 = np.array(d)[13:23]
dfreq4 = np.array(d)[14:24]
dfreq5 = np.array(d)[25:36]
dfreq6 = np.array(d)[26:37]

dsep1 = diff(dfreq2, dfreq1)
dsep2 = diff(dfreq4, dfreq3)
dsep3 = diff(dfreq6, dfreq5)
dsep4 = av(dfreq2, dfreq1)
dsep5 = av(dfreq4, dfreq3)
dsep6 = av(dfreq6, dfreq5)

ddeltfreq1 = np.array(d)[0:12]
ddeltfreq5 = np.array(d)[25:37]

ddeltsep1 = diff(ddeltfreq1, ddeltfreq5)
ddeltsep2 = av(ddeltfreq1, ddeltfreq5)

plt.figure()
plt.plot(dsep4, dsep1, "r", label="Angular degree 0")
plt.plot(dsep5, dsep2, "g", label="Angular degree 1")
plt.plot(dsep6, dsep3, "b", label="Angular degree 2")
plt.title("Log9 Large Separation as a function of Frequency")
plt.ylabel("Large Separation, $\Delta$v, $\mu$Hz")
plt.xlabel("Frequency, $\mu$Hz")
plt.legend(loc="upper left")

##########################################################################################
#star 5
efreq1 = np.array(e)[0:10]
efreq2 = np.array(e)[1:11]
efreq3 = np.array(e)[12:22]
efreq4 = np.array(e)[13:23]
efreq5 = np.array(e)[24:34]
efreq6 = np.array(e)[25:35]

esep1 = diff(efreq2, efreq1)
esep2 = diff(efreq4, efreq3)
esep3 = diff(efreq6, efreq5)
esep4 = av(efreq2, efreq1)
esep5 = av(efreq4, efreq3)
esep6 = av(efreq6, efreq5)

edeltfreq1 = np.array(e)[0:11]
edeltfreq5 = np.array(e)[24:35]

edeltsep1 = diff(edeltfreq1, edeltfreq5)
edeltsep2 = av(edeltfreq1, edeltfreq5)

plt.figure()
plt.plot(esep4, esep1, "r", label="Angular degree 0")
plt.plot(esep5, esep2, "g", label="Angular degree 1")
plt.plot(esep6, esep3, "b", label="Angular degree 2")
plt.title("Log11 Large Separation as a function of Frequency")
plt.ylabel("Large Separation, $\Delta$v, $\mu$Hz")
plt.xlabel("Frequency, $\mu$Hz")
plt.legend(loc="upper left")

#####################################################################################
# small sep

plt.figure()
plt.plot(adeltsep2, adeltsep1, "r", label= "Log3")
plt.plot(bdeltsep2, bdeltsep1, "g", label="Log5")
plt.plot(cdeltsep2, cdeltsep1, "b", label="Log7")
plt.plot(ddeltsep2, ddeltsep1, "y", label="Log9")
plt.plot(edeltsep2, edeltsep1, "k", label="Log11")
plt.title("Small Separation as a function of Frequency")
plt.ylabel(" Separation, $\delta v_{02}$, $\mu$Hz")
plt.xlabel("Frequency, $\mu$Hz")
plt.legend(loc="upper left")


average_small = []
average_small_err = []
average_large = []
average_large_err = []

average_small.append(np.average(adeltsep1))
average_small.append(np.average(bdeltsep1))
average_small.append(np.average(cdeltsep1))
average_small.append(np.average(ddeltsep1))
average_small.append(np.average(edeltsep1))

average_small_err.append(np.std(adeltsep1))
average_small_err.append(np.std(bdeltsep1))
average_small_err.append(np.std(cdeltsep1))
average_small_err.append(np.std(ddeltsep1))
average_small_err.append(np.std(edeltsep1))

lista = list(asep1)+list(asep2)+list(asep3)
listb = list(bsep1)+list(bsep2)+list(bsep3)
listc = list(csep1)+list(csep2)+list(csep3)
listd = list(dsep1)+list(dsep2)+list(dsep3)
liste = list(esep1)+list(esep2)+list(esep3)

average_large.append(np.average(lista))
average_large.append(np.average(listb))
average_large.append(np.average(listc))
average_large.append(np.average(listd))
average_large.append(np.average(liste))

average_large_err.append(np.std(lista))
average_large_err.append(np.std(listb))
average_large_err.append(np.std(listc))
average_large_err.append(np.std(listd))
average_large_err.append(np.std(liste))

stars = ["Log3", "Log5", "Log7", "Log9", "Log11"]

plt.figure()
plt.plot(mass_frac, average_small, "r.")
plt.title("Average Small Separation as a function of the \n Hydrogen Mass Fraction")
plt.errorbar(mass_frac, average_small, yerr = average_small_err, fmt = "none")
plt.ylabel("Average Small Separation, $\delta v_{02}$")
plt.xlabel("Central Hydrogen Mass Fraction, $log(X_c)$")
for i, txt in enumerate(stars):
    plt.annotate(txt, (mass_frac[i]- 0.2, average_small[i]+1.5))


plt.figure()
plt.plot(mass_frac, average_large, "r.")
plt.errorbar(mass_frac, average_large, yerr = average_large_err, fmt = "none")
plt.title("Average Large Separation as a function of the \n Central Hydrogen Mass Fraction")
plt.ylabel("Average Large Sparation, $\Delta v_{02}$")
plt.xlabel("Central Hydrogen Mass Fraction, $log(X_c$)")
for i, txt in enumerate(stars):
    plt.annotate(txt, (mass_frac[i]- 0.2, average_large[i]+1.8))

plt.figure()
plt.plot(np.average(lista), np.average(adeltsep1), "r.", label= "Log3")
plt.plot(np.average(listb), np.average(bdeltsep1), "g.", label="Log5")
plt.plot(np.average(listc), np.average(cdeltsep1), "b.", label="Log7")
plt.plot(np.average(listd), np.average(ddeltsep1), "y.", label="Log9")
plt.plot(np.average(liste), np.average(edeltsep1), "k.", label="Log11")
plt.errorbar(np.average(lista), np.average(adeltsep1), xerr = np.std(lista),yerr = np.std(adeltsep1), fmt = "none", ecolor = "r")
plt.errorbar(np.average(listb), np.average(bdeltsep1), xerr = np.std(listb),yerr = np.std(bdeltsep1), fmt = "none", ecolor = "g")
plt.errorbar(np.average(listc), np.average(cdeltsep1), xerr = np.std(listc),yerr = np.std(cdeltsep1), fmt = "none", ecolor = "b")
plt.errorbar(np.average(listd), np.average(ddeltsep1), xerr = np.std(listd),yerr = np.std(ddeltsep1), fmt = "none", ecolor = "y")
plt.errorbar(np.average(liste), np.average(edeltsep1), xerr = np.std(liste),yerr = np.std(edeltsep1), fmt = "none", ecolor = "k")
plt.title("CD Diagram")
plt.ylabel(" Separation, $\delta v_{02}$, $\mu$Hz")
plt.xlabel("Average Large Sparation, $\Delta v$")
plt.legend(loc="upper left")











