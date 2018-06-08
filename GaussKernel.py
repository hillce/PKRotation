# Script for Gaussian Kernel
# Charles Hill
# 16/05/18

import os
import matplotlib as mpl 
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.ndimage.filters import gaussian_filter1d as gd1

def GaussKernel(Path,File,sigma=20,m=1.0,b=0.0,Rescaled=False):

    if Rescaled == True:
        Cf = np.load(Path+File+"Cf_rescaled.npy")
    else:
        Cf = np.load(Path+File+"Cf.npy")

    Cf = abs(Cf)*m + b
    Cf_smooth = gd1(np.sort(abs(Cf)),sigma)

    return Cf_smooth

def Batch_GaussKernel(Path,Files,sigma=20,m=1.0,b=0.0,Rescaled=False):
    Cf_smooth = []
    for File in Files:
        print(File)
        Cf_s = GaussKernel(Path,File,sigma=sigma,m=m,b=b)
        Cf_smooth.append(Cf_s)

    return Cf_smooth

'''
m = 48193.78714638946
b = -76.26540739843597

Path = "./Contrasts/DRP1_2/"
Files = "DRP1_250_GMP-PCP_Cf.npy"
Cf_smooth = GaussKernel(Path,Files)

n, b, p = plt.hist(Cf_smooth,bins=1000,m=m,b=b)

plt.figure()
plt.plot(b[:-1],n_s)
plt.show()
'''