# Automatic Gaussian Parameters Estimation
# Charles Hill, 31/05/18

import numpy as np
from astroML.plotting import hist
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

from Fit_Gaussian import Fit_Gaussian

def Auto_Gauss(Path,File):
    Cf = np.load(Path+File+"Cf.npy")
    plt.figure
    n, bins, pat = hist(abs(Cf),bins='knuth',range=(0.0,np.max(abs(Cf))*0.5),align='left')
    plt.clf()
    Rel_Max = argrelextrema(n,np.greater,order=6)
    ctr = bins[Rel_Max]
    amp = n[Rel_Max]
    Wid = np.zeros(len(n[Rel_Max]))
    j = 0
    for idx in Rel_Max[0]:
        i = 1
        while n[idx]/2.0 < n[idx+i]:
            i += 1
        ind = idx + i
        Wid[j] = (bins[ind]-bins[idx])*2
        j += 1

    p_guess = np.zeros(len(Wid)*3)
    j = 0
    for i in range(0,len(p_guess),3):
        p_guess[i] = ctr[j]
        p_guess[i+1] = amp[j]
        p_guess[i+2] = Wid[j]
        j += 1
    
    return p_guess
