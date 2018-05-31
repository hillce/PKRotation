# Function for producing the oligomer proportion versus Concentration
# Charles Hill
# 14/05/18

import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt
import math

def Oligomer_proportion(amp):
    amp = np.array(amp)
    Olig_prop = np.zeros(amp.shape)

    for i in range(amp.shape[0]):
        Sum = 0.0
        for j in range(amp.shape[1]):
            Sum += amp[i][j]
        for k in range(amp.shape[1]):
            Olig_prop[i][k] = amp[i][k]/Sum

    return Olig_prop

def Plot_proportions(Olig_prop,Concentrations):
    Olig_prop = np.array(Olig_prop)
    no_conc = Olig_prop.shape[0]
    no_peaks = Olig_prop.shape[1]
    no_subplots = math.ceil(no_peaks/2)

    if no_subplots == 1:
        no_subplots_2 = 2
    else:
        no_subplots_2 = no_subplots 

    plt.figure()
    for i in range(no_peaks):
        plt.subplot(no_subplots,no_subplots_2,(i+1))
        peaks = []
        for j in range(no_conc):
            peaks.append(Olig_prop[j][i])
        plt.plot(Concentrations,peaks)
        plt.xlabel("Concentration")
        plt.ylabel("Oligomer Proportion")
        plt.title("Oligomer %i" % (i+1))

    plt.show()

    return


    