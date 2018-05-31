# Gaussian Fitting script
# Charles Hill
# 16/05/18

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astroML.plotting import hist
import numpy as np

def func(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)
    return y

def Fit_Gaussian(Cf,p_guess,bins,Mass=False):
    plt.figure()
    (n, bins, pathces) = hist(abs(Cf),bins=bins,normed=True,align='mid')
    print("Number of Bins:",len(bins))
    (popt, pcov) = curve_fit(func, bins[:-1], n, p0=p_guess,method='dogbox')
    if Mass == True:
        x = np.linspace(0.0,1500,1500)
    else:
        x = np.linspace(0.0,0.2,1500)
    fit = func(x, *popt)

    plt.plot(x,fit, 'b--')
    plt.yticks([])
    if Mass == True:
        plt.xlabel("Mass (kDa)",fontsize = 12, color = 'blue')
    else:
        plt.xlabel("Contrast",fontsize = 12,color='blue')
    plt.ylabel("Probability Density",fontsize = 12,color='blue')
    plt.show()

    ctr_opt = []
    amp_opt = []
    wid_opt = []

    for i in range(0, len(popt), 3):
        ctr_opt.append(popt[i])
        amp_opt.append(popt[i+1])
        wid_opt.append(popt[i+2])

    ctr_opt = np.array(ctr_opt)
    amp_opt = np.array(amp_opt)
    wid_opt = np.array(wid_opt)

    return ctr_opt, amp_opt, wid_opt, popt

def Batch_Fit(Path,Files,Rescaled=False,bins='knuth',Mass=False,m=1.0,b=0.0):
    
    ctr = []
    amp = []
    wid = []
    popt = []

    for File in Files:
        print(File)
        if Rescaled == True:
            Cf = np.load(Path+File+"Cf_rescaled.npy")
        else:
            Cf = np.load(Path+File+"Cf.npy")

        if Mass == True:
            Cf = abs(Cf)*m + b

        plt.figure()
        hist(abs(Cf),bins=bins,normed=True,align='mid')
        plt.show(block=False)

        print("Number of Gaussians:")
        No_gauss = int(input())
        print('\n')

        p_guess = []
        for i in range(No_gauss):
            print("Gaussian %i" % (i+1))
            print("Centre:")
            p_guess.append(float(input()))
            print('\n')
            print("Amplitude:")
            p_guess.append(float(input()))
            print('\n')
            print("Width:")
            p_guess.append(float(input()))
            print('\n')
        
        print("Parameter Guess:",p_guess)

        center, amplitude, width, param_opt = Fit_Gaussian(abs(Cf),p_guess,bins,Mass)

        ctr.append(center)
        amp.append(amplitude)
        wid.append(width)
        popt.append(param_opt)

    return ctr, amp, wid, popt