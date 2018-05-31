# Machine Learning of particle location

# Import packages
import os

import iscat.read as isc_rd
import iscat.cut as isc_cut
import iscat.stats as isc_ss
import iscat.detect as isc_dt
import iscat.thumbnail as isc_tl
import iscat.peak_quality as isc_py
import iscat.contrast as isc_ct
import iscat.reflectivity as isc_ry
import iscat.peak as isc_pk
import iscat.psf as isc_pf
import iscat.fit as isc_ft
import iscat.utils as isc_us
import iscat.stack as isc_sk
import iscat.plot as isc_pt

import cv2
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from Multi_slice import multi_slice_viewer

# Color scale settings for imshow of peaks
#amp = 1E-2
amp = 4E-3
kwargs = {"vmin": -amp, "vmax": amp, "cmap" : "RdBu_r"}

# Specify TDMS file to be analysed
DIR_SAMPLE = r"C:/Users/PKGroup/Documents/CH/Charlie/2018/May/"
FILE_SAMPLE = r"03/ADH10_0/cam1/event2.tdms"

filename_sample = os.path.join(os.path.expanduser(DIR_SAMPLE),FILE_SAMPLE)
assert os.path.exists(filename_sample)

DETECTOR_FULL_WELL_CAPACITY = 32513
DETECTOR_BINNING_TIME       = 10 # Generally set to 5, remember to take note during experiment
DETECTOR_BINNING_PIXEL      = 5**2 # Generally set to 5**2, as above
DETECTOR_GREY_LEVELS_DATA   = 4096

# Read data from file
frames_raw = isc_rd.read_data(filename_sample, verbose=False)

# Rescale raw ADU to photoelectron counts
frames = isc_rd.rescale_to_electrons(frames_raw, 
                                     full_well_capacity=DETECTOR_FULL_WELL_CAPACITY, 
                                     binning_time=DETECTOR_BINNING_TIME, 
                                     binning_pixel=DETECTOR_BINNING_PIXEL,
                                     grey_levels_data=DETECTOR_GREY_LEVELS_DATA)

print("Read %i frames of shape %i x %i" % (frames.shape[0], frames.shape[1], frames.shape[2]))

multi_slice_viewer(frames)
plt.show()

iscat_img = np.zeros(frames.shape)
for i in range(10,frames.shape[0]-10):
    temp_imgs = np.zeros((20,frames.shape[1],frames.shape[2]))
    for j in range(temp_imgs.shape[0]):
        temp_imgs[j][:][:] = frames[i-j+10][:][:]
    med_img = np.median(temp_imgs,axis=0)
    iscat_img[i][:][:] = frames[i][:][:]/med_img - 1

multi_slice_viewer(iscat_img)
plt.show()

iscat_img_avg = np.zeros((frames.shape[0]/10,frames.shape[1],frames.shape[2]))
for i in range(10,frames.shape[0],10):
    temp_imgs = np.zeros((11,frames.shape[1],frames.shape[2]))
    for j in range(temp_imgs.shape[0]):
        temp_imgs[j][:][:] = frames[i-j+5][:][:]
