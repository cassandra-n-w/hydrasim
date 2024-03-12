#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 12:39:30 2024

@author: cassie
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col


file_water = "oh2o1600.fits"
file_nowater = "oh2o1600_nowater.fits"
noiselev = 3e-3

imsize_pixels = 1000

# do plotting of radially averaged observation
beamsize = 2.4e-1 # beamsize in arcseconds, obtained by averaging the major and minor axes

spacing = 10
means = []
radii = np.arange (0.01,imsize_pixels/2, spacing)

radtodeg = 180/np.pi
degtoarcsec = 3600

#pixsize = abs(headout['incr'][0]) * radtodeg * degtoarcsec
pixsize = 0.0049921681
radii_arcsec = radii*pixsize
radii_beams = radii_arcsec/beamsize

numseeds = 15

means = np.load("means0.npy")
multimean_water = np.zeros((len(means), numseeds))
multimean_nowater = np.zeros((len(means), numseeds))

for seedmod in np.arange(0,numseeds):
    
    means = np.load(str(int(noiselev*1000)) + "mjy/means" + str(seedmod) + file_water + ".npy")
    multimean_water[:,seedmod] = means
    means = np.load(str(int(noiselev*1000)) + "mjy/means" + str(seedmod) + file_nowater + ".npy")
    multimean_nowater[:,seedmod] = means


mean_water_stddev = np.std(multimean_water, axis=1)
mean_water_avg = np.mean(multimean_water, axis=1)

mean_nowater_stddev = np.std(multimean_nowater, axis=1)
mean_nowater_avg = np.mean(multimean_nowater, axis=1)

for seedmod in np.arange(0,numseeds):
    color_water = col.hsv_to_rgb(((204 -15+ 2*seedmod)/360, 83/100, 70.5/100))
    color_nowater = col.hsv_to_rgb(((28 -15+ 2*seedmod)/360, 97/100, 100/100))
    means = multimean_water[:,seedmod]
    plt.errorbar(radii_arcsec, means, ls="--", yerr=mean_water_stddev, color=color_water)
    means = multimean_nowater[:,seedmod]
    plt.errorbar(radii_arcsec, means, yerr=mean_nowater_stddev, color=color_nowater)
    
range_min = range(8,20)
range_max = range(15, 30)
range_of_interest = range(11,25)

#signal_noiseless = np.max(means_noiseless[range_max]) - np.min(means_noiseless[range_min ])
# signal_noisy = np.max(mean_avg[range_max ]) - np.min(mean_avg[range_min ])

# average_noise = np.mean(mean_stddev[range_of_interest])


# #print("noiseless signal " + str(signal_noiseless*1000) + "mJy/beam")
# print("noisy average signal " + str(signal_noisy*1000) + "mJy/beam")

# print("noise in region of interest " + str(average_noise*1000) + "mJy/beam")
# print("SNR " + str(10*np.log10(signal_noisy/average_noise)) + "dB")
# print("SNR " + str(signal_noisy/average_noise) + " (linear)")

plt.title("Azimuthally Averaged Intensity Comparison of Water vs Continuum, Noise = " + str(noiselev*1000) + "mJy/beam")
plt.xlabel("Distance from disk center (arcsec)")
plt.ylabel("Average Intensity (Jy/beam)")
plt.xlim(0,max(radii*pixsize))
plt.ylim(-0.01, 0.03)
plt.legend(["With Water", "Without water (continuum only)"])
plt.show()


plt.errorbar(radii_arcsec, mean_water_avg, mean_water_stddev)
plt.errorbar(radii_arcsec, mean_nowater_avg, mean_nowater_stddev)
plt.title("Azimuthally Averaged Intensity Comparison of Water vs Continuum, Noise = " + str(noiselev*1000) + "mJy/beam")
plt.xlabel("Distance from disk center (arcsec)")
plt.ylabel("Average Intensity (Jy/beam)")
plt.xlim(0,max(radii*pixsize))
plt.ylim(-0.01, 0.03)
plt.legend(["With Water", "Without water (continuum only)"])
plt.show()



# plt.plot(radii_arcsec, mean_stddev/errors)
# plt.title("Ratio of Monte Carlo Noise to Statistically Estimated Noise")
# plt.xlabel("Distance from disk center (arcsec)")
# plt.ylabel("Noise Ratio (monte carlo estimate / statistical estimate)")
# plt.show()

