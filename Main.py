# Main
# Charles Hill, 31/05/18

import os
import csv
import scipy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astroML.density_estimation import knuth_bin_width, freedman_bin_width, scotts_bin_width
from astroML.plotting import hist

# Own Functions and files
from Poly_Fig import Polygon_figure, Gauss_sum
from Oligomer_proportions import Oligomer_proportion, Plot_proportions
from Fit_Gaussian import Batch_Fit
from Calibration import Calibration_Line, Convert_mass, Convert_fits
from Normalisation import normalisation_max, normalisation_max_1

#Nick's Fit program
from EMfit import fit

Buffer = "DRP1/"
Files = ["ADH_STD_","BetaA_STD_"]
m,b = Calibration_Line(Buffer,Files)

Path = "C:/Users/PKGroup/iscat/iscat_jupyter/Contrasts/DRP1_2/"
Files = ["DRP1_250_","DRP1_250_GMP-PCP_","DRP1_250_GDP_","DRP1_500_"]
handles = ["250nM","250nM with GMP-PCP","250nM with GDP","500nM"]
#Files = ["DRP1_250_GMP-PCP_","DRP1_250_GDP_"]
#handles = ["250nM with GMP-PCP","250nM with GDP"]
ctr, amp, wid, popt = Batch_Fit(Path,Files,Mass=True,m=m,b=b)

print("Optimised Parameters:",popt)

x, y = normalisation_max_1(popt,4,1500,x_point=1000)

plt.figure()
for i in range(y.shape[0]):
    plt.plot(x,y[i,:])
plt.legend(handles)
plt.xlabel("Mass/kDa")
plt.ylabel("P(M)")
plt.show()

#Conc = [10,50,100,500,1000]
#Polygon_figure(Conc,popt)

print("Hello")


print("Bye")