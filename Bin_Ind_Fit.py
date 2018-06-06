# iSCAMS class file.
# Charles Hill
# 03/06/18

import numpy as np
import math
import matplotlib.pyplot as plt 
from astroML.plotting import hist
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema
from scipy.ndimage.filters import gaussian_filter1d as gd1
from Poly_Fig import Gauss_sum
from iSCAMS_class import iSCAMS

from scipy.odr import ODR, Model, Data, RealData
from astroML.plotting import hist

Path = "C:/Users/Charles Hill/OneDrive - Nexus365/Contrasts/DRP1_2/"
File = "DRP1_250_GMP-PCP_"
Cf = np.sort(abs(np.load(Path+File+"Cf.npy")))
'''
for i in range(0,20,4):
    print(i)
    Cf_temp = Cf*10**i
    DRP1 = iSCAMS(Cf_temp,Conc="250nM",Protein="DRP1",Buffer="HEPES",Nucleotide="GMP-PCP",Mass=False,c_range=(0.0,np.max(Cf_temp)))
    DRP1.Auto_Gauss()
    DRP1.Fit_Gaussian()

def func(beta, x): # Multiple Gaussian function for scipy.optimize.curve_fit
    y = np.zeros_like(x)
    for i in range(0, len(beta), 3):
        ctr = beta[i]
        amp = beta[i+1]
        wid = beta[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)
    return y
'''
def func(x, *params): # Multiple Gaussian function for scipy.optimize.curve_fit
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)
    return y

guess = [0.0055,4000,0.001,0.0105,700,0.002,0.0145,200,0.002,0.019,60,0.002,0.025,15,0.002]
LS_output = np.zeros(len(guess))
ODR_output = np.zeros(len(guess))

X_tot = np.zeros((5,500))
Y_tot = np.zeros((5,500))
plt.figure()
for i in range(5):
    Cf_temp = Cf*10**i
    n, bins, p = hist(Cf_temp,bins='knuth')
    #data = Data(bins[:-1],n)
    #model = Model(func)

    new_guess = []
    for j in range(0,len(guess),3):
        new_guess.append(guess[j]*10**i)
        new_guess.append(guess[j+1])
        new_guess.append(guess[j+2]*10**i)

    print(new_guess)

    popt, pcov = curve_fit(func,bins[:-1],n,p0=new_guess)
    xn = np.linspace(0,0.03*(10**i),500)
    yn = func(xn,*popt)

    X_tot[i,:] = xn/(10**i)
    Y_tot[i,:] = yn

    plt.plot(xn,yn,label='Contrast to the power of %i'% (i))
    '''
    odr = ODR(data, model, new_guess)
    odr.set_job(fit_type=2)
    output = odr.run()
    
    xn = np.linspace(0,0.03*(10**i),500)
    yn = func(output.beta,xn)
    plt.plot(xn,yn,'k-',label='least sq')
    for j in range(0,len(guess),3):
        LS_output[j]

    odr.set_job(fit_type=0)
    output = odr.run()
    yn = func(output.beta,xn)
    plt.plot(xn,yn,'g-',label='odr')
    ODR_output += output.beta
    '''
plt.xlabel("Contrast")
plt.legend(loc=0)
plt.show()

plt.figure()
for i in range(5):
    plt.plot(X_tot[i,:],Y_tot[i,:],label='Contrast to the power of %i'% (i))
plt.xlabel("Contrast")
plt.legend(loc=0)
plt.show()




'''
x = np.linspace(1,len(Cf),len(Cf))

Cf_diff = np.diff(np.diff(Cf))

plt.figure()
plt.subplot(1,2,1)
plt.plot(x,Cf)
plt.subplot(1,2,2)
plt.plot(x[:-2],Cf_diff)
plt.show()

mask = (-0.001 <= Cf_diff) & (Cf_diff <= 0.001)
mask = np.append(mask,(False,False))
print(mask)
New_Cf = Cf[mask]

plt.figure()
plt.subplot(1,2,1)
plt.hist(Cf,bins=1000)
plt.subplot(1,2,2)
plt.hist(New_Cf,bins=1000)
plt.show()
'''