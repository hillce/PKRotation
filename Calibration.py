# Calibration python script
# Charles Hill
# 16/05/2018

from Fit_Gaussian import Batch_Fit
import numpy as np
from numpy.polynomial.polynomial import polyfit
import matplotlib.pyplot as plt

def Calibration_Line(Path,Buffer,Files):
    Path = Path
    Buffer = Buffer
    Files = Files

    STD_DICT = {"ADH_STD_":(77.8,152.3),"BetaA_STD_":(60.7),"ProA_STD_":(44),"BSA_STD_":(66.6)}

    ctr, amp, wid, popt = Batch_Fit(Path+Buffer,Files)

    print(ctr)

    Mass = []
    Centres = []
   
    i = 0
    for Protein in Files:
        Peaks = STD_DICT[Protein]
        if type(Peaks) == tuple:
            for k in range(len(Peaks)):
                Mass.append(Peaks[k])
        else:
            Mass.append(Peaks)
        for j in range(len(ctr[i])):
            Centres.append(ctr[i][j])
        i += 1

    b, m = polyfit(Centres,Mass, 1)

    x = np.arange(0.001,0.01,0.001)
    y = b + m*x

    plt.plot(Centres, Mass, '.')
    plt.plot(x, y, '-')
    plt.show()

    return m, b

def Convert_mass(Path,Files,m,b,Rescaled=False):

    Cf_mass_dict = {}
    for File in Files:
        print(File)
        if Rescaled == True:
            Cf = np.load(Path+File+"Cf_rescaled.npy")
        else:
            Cf = np.load(Path+File+"Cf.npy")

        Cf_mass_dict[File] = abs(Cf)*m + b
    
    return Cf_mass_dict

def Convert_fits(m,b,ctr,amp,wid,popt):
    
    ctr_mass = ctr*m + b
    amp_mass = amp*m + b
    wid_mass = wid*m + b
    popt_mass = popt*m + b

    return ctr_mass, amp_mass, wid_mass, popt_mass