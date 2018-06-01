# iSCAMS class file.
# Charles Hill
# 01/06/18
import numpy as np
import matplotlib.pyplot as plt 
from astroML.plotting import hist
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema
from scipy.ndimage.filters import gaussian_filter1d as gd1

###############################################################################################################################
#### The iSCAMS class contains all the functions and info that you would want to extract from the contrast of the fitted   ####
#### peaks from Max's Program (see C:\Users\Charles Hill\iscat\iscat_jupyter\3_iscat_pipeline.ipynb). Not to be used with  ####
#### multiple contrasts. There is a separate set of functions for dealing with the outputs from this class collectively    ####
###############################################################################################################################

class iSCAMS:
    def __init__(self,Cf,Conc='',Protein = '',Buffer = '',Nucleotide = '',c_range=(0.0,0.15),m_range=(0.0,1500),Mass=True,bin_type='knuth',order=3):
        self.contrast = abs(Cf)
        self.instances = len(Cf)
        self.mass_data = []
        self.popt = []
        self.p_guess = []
        self.x = []
        self.fit = []
        self.bins = bin_type
        self.Mass = Mass
        self.c_range = c_range
        self.m_range = m_range
        self.order = order
        self.contrast_smooth = []
        self.mass_smooth = []
        self.Conc = Conc
        self.Protein = Protein
        self.Buffer = Buffer
        self.Nucleotide = Nucleotide

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
        if self.Mass == True:
            n, bins, pathces = hist(self.mass_data,bins=self.bins,normed=True,align='mid')
        else:
            n, bins, pathces = hist(self.contrast,bins=self.bins,normed=True,align='mid')
        while True:
            try:
                self.popt, pcov = curve_fit(self.func, bins[:-1], n, p0=self.p_guess)
                break
            except RuntimeError:
                # If the curvefit fails, then instead of returning the errors, it recalls the Auto_gauss with a different order
                print("Error - curve_fit failed, rerunning Auto_Gauss")
                self.order += 1
                self.Auto_Gauss()
                

        if self.Mass == True:
            self.x = np.linspace(self.m_range[0],self.m_range[1],1500)
        else:
            self.x = np.linspace(self.c_range[0],self.c_range[1],1500)
        
        self.fit = self.func(self.x, *self.popt)

        plt.plot(self.x,self.fit, 'b--')
        plt.yticks([])
        if self.Mass == True:
            plt.xlabel("Mass (kDa)",fontsize = 12, color = 'blue')
        else:
            plt.xlabel("Contrast",fontsize = 12,color='blue')
        plt.ylabel("Probability Density",fontsize = 12,color='blue')
        plt.show()

    def Auto_Gauss(self):
        plt.figure
        if self.Mass == True:
            n, bins, pat = hist(self.mass_data,bins=self.bins,align='left')
        else:
            n, bins, pat = hist(self.contrast,bins=self.bins,align='left')

        Rel_Max = argrelextrema(n,np.greater,order=self.order)
        ctr = bins[Rel_Max]
        amp = n[Rel_Max]

        temp_idx = []
        for i in range(len(amp)):
            if amp[i] <= 3.0:
                temp_idx.append(i)
        amp = np.delete(amp,temp_idx)
        ctr = np.delete(ctr,temp_idx)
        Rel_Max = np.delete(Rel_Max,temp_idx)

        Wid = np.zeros(len(amp))
        j = 0
        print(Rel_Max)
        for idx in Rel_Max:
            i = 1
            while n[idx]/2.0 < n[idx+i]:
                i += 1
            ind = idx + i
            Wid[j] = (bins[ind]-bins[idx])*2
            j += 1

        self.p_guess = np.zeros(len(Wid)*3)
        j = 0
        for i in range(0,len(self.p_guess),3):
            self.p_guess[i] = ctr[j]
            self.p_guess[i+1] = amp[j]
            self.p_guess[i+2] = Wid[j]
            j += 1

        plt.plot(ctr,amp,'r.')
        plt.plot(ctr-Wid/2.0,amp/2.0,'g.')
        plt.plot(ctr+Wid/2.0,amp/2.0,'g.')
        plt.show()

    def Manual_Gauss(self):
        plt.figure
        if self.Mass == True:
            hist(self.mass_data,bins=self.bins,align='left')
        else:
            hist(self.contrast,bins=self.bins,align='left')
        plt.show()

        print("Number of Gaussians:")
        No_gauss = int(input())
        for i in range(No_gauss):
            print("Gaussian %i" % (i+1))
            print("Centre:")
            self.p_guess.append(float(input()))
            print("Amplitude:")
            self.p_guess.append(float(input()))
            print("Width:")
            self.p_guess.append(float(input()))

    def GaussKernel(self,sigma=20):
        if self.Mass == True:
            self.mass_smooth = gd1(np.sort(self.mass_data),sigma)
        else:
            self.contrast_smooth = gd1(np.sort(self.contrast),sigma)

    def Plot_hist(self):
        if self.Mass == True:
            plt.figure()
            hist(self.mass_data,bins=self.bins,align='left',label=self.Protein+"("+self.Conc+","+self.Buffer+","+self.Nucleotide+")")
            plt.legend()
            plt.xlabel("Mass (kDa)")
            plt.ylabel("Particle Frequency")
            plt.show()
        else:
            plt.figure()
            hist(self.contrast,bins=self.bins,align='left',label=self.Protein+"("+self.Conc+","+self.Buffer+","+self.Nucleotide+")")
            plt.legend()
            plt.xlabel("Contrast")
            plt.ylabel("Particle Frequency")
            plt.show()

    


    
