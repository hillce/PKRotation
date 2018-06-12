# Calibration of setups:
# Charles Hill
# 08/06/18

import numpy as np
from numpy.polynomial.polynomial import polyfit
import matplotlib.pyplot as plt
from iSCAMS_class import iSCAMS
from astroML.plotting import hist
from scipy.optimize import curve_fit

class iSCAMS_calibrate():
    def __init__(self,Cf,Protein,bins="knuth"):
        self.STD_DICT= {"ADH":(77.8,152.3),"BetaA":(60.7),"ProA":(44),"BSA":(66.6)}
        self.Protein = Protein
        self.contrast = abs(Cf)
        self.bins = bins
        self.popt = []
        self.p_guess = []
        self.ctr = []

    def func(self, x, *params): # Multiple Gaussian function for scipy.optimize.curve_fit
        y = np.zeros_like(x)
        for i in range(0, len(params), 3):
            ctr = params[i]
            amp = params[i+1]
            wid = params[i+2]
            y = y + amp * np.exp( -((x - ctr)/wid)**2)
        return y

    def Fit_Gaussian(self): # Fits multiple Gaussians according to the parameter guess input (either manual or auto)
        plt.figure()
        n, bins, pathces = hist(self.contrast,bins=self.bins,normed=True,align='left')
        self.popt, pcov = curve_fit(self.func, bins[:-1], n, p0=self.p_guess)

        for i in range(0, len(self.popt), 3):
            self.ctr = self.popt[i]

    def Manual_Gauss(self): # Seeds the parameters, based on manual input.
        plt.figure
        hist(self.contrast,bins=self.bins,align='left')
        plt.show()

        print("Number of Gaussians:")
        No_gauss = int(input())
        for i in range(No_gauss):
            print("Gaussian %i" % (i+1))
            print("Centre:")
            np.append(self.p_guess,float(input()))
            print("Amplitude:")
            np.append(self.p_guess,float(input()))
            print("Width:")
            np.append(self.p_guess,float(input()))