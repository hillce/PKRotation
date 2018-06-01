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

DRP1 = iSCAMS(Cf,Mass=False)

print(DRP1.c_range)
print(DRP1.contrast)

print(DRP1.p_guess)

DRP1.Auto_Gauss()

print(DRP1.p_guess)