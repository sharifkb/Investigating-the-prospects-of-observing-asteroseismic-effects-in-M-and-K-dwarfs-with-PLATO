import numpy as np
import matplotlib.pyplot as plt
import  astropy.constants as c
import matplotlib as mpl
mpl.rc('font', family='serif')

#constants
R_sun = c.R_sun.to("km").value 
T_eff_sun = 5777
g_sun = 274
tau_sun = 200e-6
V_max_sun = 3100
delta_v_sun = 134.8
A_sun = 2.1
D_sun = 1.493
epsilon=1.5
sum_vis = 0.5 + 1.5 + 0.04

#frequency space
#n = np.linspace(18,21,4)
n = np.linspace(1,100,100)
x = np.linspace(1,15000, 1000000)


def spectrum (T_eff, g, M, L, vmag, log = True):     
        
    #granulation    
    R = np.sqrt(L*c.L_sun.value/(4*np.pi*c.sigma_sb.value*T_eff**4))*10**(-3)    
    tau_c = tau_sun*(g_sun/g)*(T_eff/T_eff_sun)**(1/2)    
    sigma_c = 50000/(c.R_sun.value/(1000*10**3))*(T_eff/T_eff_sun)*R/(M*R_sun)  
    gran = 4*np.sqrt(2)*tau_c*sigma_c**2/(1+(2*np.pi*x*tau_c)**4)  
    
    #finding the mode locations
    delta_v = delta_v_sun*M**(1/2)*(R/R_sun)**(-3/2)
    sigma = 2.123*delta_v
    D = D_sun*(delta_v/delta_v_sun)   
    
    modes0=[]
    modes1=[]
    modes2=[]
    modes3=[]
    for i in range(len(n)):
        modes0.append(delta_v*(n[i]    +epsilon))
        modes1.append(delta_v*(n[i]+1/2+epsilon)-D*2)  
        modes2.append(delta_v*(n[i]+  1+epsilon)-D*6)  
        modes3.append(delta_v*(n[i]+3/2+epsilon)-D*12)   
    
    #introducing the linewidth and mode heights
    linewidth = 0.07+0.91*(T_eff/5777)**15.3
    beta = 1 - np.exp(-(8907*L**-0.093 - T_eff)/1250)
    V_max = V_max_sun*(g/(g_sun))*(T_eff/T_eff_sun)**(-1/2)    
    A_max = beta*A_sun*L/M*(T_eff/T_eff_sun)**(-2)      
    Gauss = ((A_max**2*sum_vis)/(2*delta_v))*np.exp(-0.5*((x-V_max)/sigma)**2)          
    height_new = (Gauss*2*delta_v)/(linewidth*np.pi*sum_vis)
    
    a0 = np.zeros(len(x))
    a1 = np.zeros(len(x))
    a2 = np.zeros(len(x))
    a3 = np.zeros(len(x))
    for i in range(0,len(n)):
        a0+=height_new     /(1+((x - np.array(modes0[i]))/(linewidth/2))**2)
        a1+=height_new*1.5 /(1+((x - np.array(modes1[i]))/(linewidth/2))**2)   
        a2+=height_new*0.5 /(1+((x - np.array(modes2[i]))/(linewidth/2))**2)   
        a3+=height_new*0.04/(1+((x - np.array(modes3[i]))/(linewidth/2))**2)           
    
    #adding shotnoise
    vmag_value = np.full((len(x)),vmag)
    shotnoise = 18.0*10**(-0.4*(11.0-vmag_value))

    sig0 = a0 + gran + shotnoise
    sig1 = a1 + gran + shotnoise
    sig2 = a2 + gran + shotnoise
    sig3 = a3 + gran + shotnoise

    fig1 = plt.figure(figsize=(8,5.5))
    if log:
        plt.plot(x, np.log10(sig0), label = "l = 0", color = "red")
        plt.plot(x, np.log10(sig1), label = "l = 1", color = "green")
        plt.plot(x, np.log10(sig2), label = "l = 2", color = "blue")
        plt.plot(x, np.log10(sig3), label = "l = 3", color = "black")
    else:
        plt.plot(x, sig0, label = "l = 0", color = "red")
        plt.plot(x, sig1, label = "l = 1", color = "green")
        plt.plot(x, sig2, label = "l = 2", color = "blue")
        plt.plot(x, sig3, label = "l = 3", color = "black")       
    plt.ylabel(r"Power ($ppm^2/ \mu Hz$)", fontsize = 20)
    plt.xlabel(r"Frequency ($\mu Hz$)", fontsize = 20)   
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.legend(fontsize = 24)
    plt.tight_layout()
    plt.show()
   
    return(fig1)

###### M DWARFS
catalogue = np.load(r"C:\Users\shari\OneDrive\Documents\Work\3rd Yr\Group studies\Code\Data array creating code\Mdwarf\M_dwarves_plat_1_final.npy")

###### K DWARFS
#catalogue = np.load(r"C:\Users\shari\OneDrive\Documents\Work\3rd Yr\Group studies\Code\Data array creating code\Kdwarf\K_dwarves_plat_1_final.npy")

logg = catalogue["logg"]
logTe = catalogue["logTe"]
logL = catalogue["logL"]
Mass = catalogue["Mass"]
vmag = catalogue["V"]

#CHANGE I FROM "1" TO CHOOSE STARS TO BE PLOTTED
i= [0]   
#spectrum(10**logTe[i], 10**logg[i]/100, Mass[i], 10**logL[i], vmag[i])

#for the sun:
spectrum(T_eff_sun, g_sun, 1, 1, 4.83)