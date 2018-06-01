# iSCAMS Main
# Charles Hill
# 01/06/18

import os
import csv
import scipy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from iSCAMS_class import iSCAMS

Path = "C:/Users/Charles Hill/OneDrive - Nexus365/Contrasts/DRP1_2/"
File = "DRP1_250_GMP-PCP_"
Cf = np.load(Path+File+"Cf.npy")

print(Cf)

DRP1 = iSCAMS(Cf,Conc="250nM",Protein="DRP1",Buffer="HEPES",Nucleotide="GMP-PCP",Mass=False)
DRP1.Auto_Gauss()
DRP1.Fit_Gaussian()

'''
plt.figure()
for i in range(1,10):
    plt.subplot(3,3,i)
    DRP1.GaussKernel(sigma=i*5)
    n,b,p = plt.hist(DRP1.contrast_smooth,bins=1000)
    plt.plot(b[:-1],n)
plt.show()
'''
DRP1.Plot_hist()