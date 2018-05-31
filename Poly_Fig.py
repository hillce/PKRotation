# Polygon Figure plotting function
# Charles Hill
# 16/05/18

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection

def Gauss_sum(x, params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)
    return y

def Polygon_figure(Conc,popt,step=0.0001,Con_lim=(0,0.01),G_max=600):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    idx = np.arange(0,len(Conc),1)
    popt = np.array(popt)
    verts = []

    for z, i in zip(Conc,idx):
        x = np.arange(Con_lim[0],Con_lim[1],step)
        y = Gauss_sum(x,popt[i][:])
        verts.append(list(zip(x,y))) # generates a list of (x,y) tuples for each Concnetration

    poly = PolyCollection(verts)
    poly.set_alpha(0.7) # Sets the trasnparency of the polygon plots
    ax.add_collection3d(poly, zs=Conc, zdir='y')

    ax.set_xlabel('Contrast')
    ax.set_xlim3d(Con_lim[0],Con_lim[1])
    ax.set_ylabel('Concentration')
    ax.set_ylim3d(Conc[0]-10,Conc[-1]+10)
    ax.set_zlabel('Gaussian')
    ax.set_zlim3d(0, G_max)

    plt.show()

    return