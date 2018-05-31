# Normalisation Function Set - Takes fitted parameters 
# Charles Hill, 29/05/2018

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Own Functions and files
from Poly_Fig import Gauss_sum

#popt = [[4.46332209e-03, 3.59607523e+02, 1.09600909e-03, 8.67526548e-03,
#        6.13666314e+01, 2.15139609e-03, 1.39021353e-02, 1.28183094e+01,
#        1.54320236e-03],[4.94425260e-03, 3.89300040e+02, 1.06471322e-03,
#        9.85291612e-03, 6.46673908e+01, 1.35599956e-03, 1.49368160e-02,
#        1.85406260e+01,1.32471748e-03, 1.99515624e-02, 5.32368814e+00, 1.59006276e-03]]

def normalisation_max(popt,num_files,range_max):
    x = np.linspace(0,range_max,num=100)
    store = np.zeros((num_files,len(x)))

    for i in range(num_files):
        params = popt[i]
        y = Gauss_sum(x,params)
        store[i,:] = y

    plt.figure()
    for i in range(num_files):
        plt.plot(x,store[i,:])
    plt.show()

    max_loc = np.argmax(store)
    col = max_loc//store.shape[1]
    row = max_loc % store.shape[1]
    
    for i in range(num_files):
        local_max = np.argmax(store[i,:])
        ratio = store[col,row]/store[i,local_max]
        store[i,:] = store[i,:]*ratio

    plt.figure()
    for i in range(num_files):
        plt.plot(x,store[i,:])
    plt.show()

    return x, store

def normalisation_max_1(popt,num_files,range_max,x_point=100):
    x = np.linspace(0,range_max,num=x_point)
    store = np.zeros((num_files,len(x)))

    for i in range(num_files):
        params = popt[i]
        y = Gauss_sum(x,params)
        store[i,:] = y

    plt.figure()
    for i in range(num_files):
        plt.plot(x,store[i,:])
    plt.show()
    
    for i in range(num_files):
        local_max = np.argmax(store[i,:])
        ratio = 1.0/store[i,local_max]
        store[i,:] = store[i,:]*ratio

    plt.figure()
    for i in range(num_files):
        plt.plot(x,store[i,:])
    plt.show()

    return x, store