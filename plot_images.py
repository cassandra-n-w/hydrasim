# This file is part of Hydrasim, a small piece of software intended to 
# simulate observations of TW Hydra.

# Hydrasim is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# Hydrasim is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with Hydrasim. If not, see <https://www.gnu.org/licenses/>. 


from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import math
import os

from radmc3dPy.image import *    
from radmc3dPy.analyze import *  
from radmc3dPy.natconst import * 
from radmc3dPy import *

#
# First set up the model with
#
#   python problem_setup.py
#
# Then make sure to have run
#
#   radmc3d mctherm
#
#os.system('radmc3d mctherm')

#
# Now make the image
#
#   radmc3d image lambda 10 incl 60
#

#we use inclination = 7 degrees, cite
# 



c = 1e6 * 299792458; #um/s
spread = 2 * 1e3 * 1e6; #10km/s, in um/s
df = spread/c

molecule_name = "oh2o"

if molecule_name == "co":
    co = readMol('co')
    # note: CO 3-2 line is index 3, 345.796 GHz
    lin_idx = 3;
    mol_idx = 1;
    freq = 345.7959899e9;
    print("Setting co as primary molecule")
    #freq += 20;

# ortho water ground transition, line 1
if molecule_name == "oh2o":
    h2o = readMol('oh2o')
    lin_idx = 1;
    mol_idx = 1;
    freq = 556.93607e9;
    print("Setting oh2o as primary molecule")

# # para water ground transition, line 1
if molecule_name == "ph2o":
    h2o = readMol('ph2o')
    lin_idx = 1;
    mol_idx = 1;
    freq = 1113.34306e9;
    print("Setting ph2o as primary molecule")

lamb = c/freq ;# in microns
#lamb=10
pixel_count = 1600
cnlam = 15;
lambmin = lamb * (1-df)
lambmax = lamb * (1+df)

#TW hya information
twhya_dist_pc = 60.1 #dist from tw hya to earth in parsecs
twhya_coord = '11h01m51.9054s -34d42m17.0316s' # coordinates of TW hya
twhya_incl = 7


clam = np.linspace(lambmin, lambmax, cnlam)

with open('camera_wavelength_micron.inp','w+') as f:
    f.write('%d\n'%(cnlam))
    for value in clam:
        f.write('%13.6e\n'%(value))
        
data = analyze.readData(ddens=True, gdens=True, gtemp=True, dtemp=True)    

#%% write out to fits file: step1
out_widthkms = 2
out_nlam = 31
makeImage(npix = pixel_count, 
          incl = twhya_incl, 
          sizeau=480, 
          #lambdarange=(lambmin, lambmax), 
          linenlam=out_nlam,
          widthkms = out_widthkms,
          iline = lin_idx)
#%% write to fits file: step 2
out_bandwidth = (freq/1e6) * out_widthkms/(c/1e9) / (out_nlam)
im = readImage()  
#plotImage(im,log=True,maxlog=6,au=True,bunit='inu',cmap=cm.hot,ifreq=0)
im.writeFits(fname=molecule_name+str(pixel_count)+".fits", 
             dpc = twhya_dist_pc, 
             coord = twhya_coord, 
             bandwidthmhz= out_bandwidth,
             casa=False,
             ifreq=0)
#plotImage(im,log=True,maxlog=6,au=True,bunit='inu',cmap=cm.hot,ifreq=1, arcsec = True, dpc = twhya_dist_pc)
#%% Plot optical density

data.getTau(wav=lamb*0.99)
taux = data.taux
#plt.plot(data.taux)
#c = plt.contour(data.grid.x/natconst.au, np.pi/2.-data.grid.y, data.taux[:,:,0].T, [1.0],  colors='b', linestyles='solid')
#plt.clabel(c, inline=1, fontsize=10)

#%%
#os.system('radmc3d image lambda 10 incl 7 sizeau 100')
#os.system('radmc3d image lambda '  + str(lamb) + ' incl 7 sizeau 300')
makeImage( incl = 90, phi = 0, wav = lamb, sizeau = 400, npix = pixel_count)

#
# to get the image
#
im   = readImage()
dum  = im.image.shape
nx   = dum[0]
nx2  = nx//2
fig = plt.figure()
fig.set_size_inches(14,10)
plotImage(im,log=True,maxlog=6,au=True,bunit='inu',cmap=cm.hot)
plt.show()

#%%

#os.system('radmc3d image lambda 10 incl 7 sizeau 100')
#os.system('radmc3d image lambda '  + str(lamb) + ' incl 7 sizeau 300')
makeImage( incl = 0, phi = 0, wav = lamb, sizeau = 550, npix = pixel_count)

#
# to get the image
#
im   = readImage()
dum  = im.image.shape
nx   = dum[0]
nx2  = nx//2
fig = plt.figure()
fig.set_size_inches(14,10)
plotImage(im,log=True,maxlog=6,au=True,bunit='inu',cmap=cm.hot)
plt.show()

#%%
#plot spectrum incl 90
#makeImage( incl = 7, phi = 0, lambdarange = (lamb * (1-df), lamb * (1+df)), nlam = 10, sizeau = 300, npix = pixel_count)
#command = 'radmc3d spectrum incl 7 sizeau 300 lambdarange ' + str(lambmin)  + ' ' + str(lambmax) + ' nlam ' + str(fcount)
command = 'radmc3d spectrum incl 90 sizeau 1000 iline ' + str(lin_idx) + ' imolspec ' + str(mol_idx) + ' widthkms 5. linenlam 100'#'lambdarange 5. '#'nlam 100'
print(command)
os.system(command)


plotSpectrum(readSpectrum(), mol=co, ilin = lin_idx, fnu = True, dpc = 60.1)



#%%
#plot spectrum incl 0
#makeImage( incl = 7, phi = 0, lambdarange = (lamb * (1-df), lamb * (1+df)), nlam = 10, sizeau = 300, npix = pixel_count)
#command = 'radmc3d spectrum incl 7 sizeau 300 lambdarange ' + str(lambmin)  + ' ' + str(lambmax) + ' nlam ' + str(fcount)
command = 'radmc3d spectrum incl 5 sizeau 500 iline ' + str(lin_idx) + ' imolspec ' + str(mol_idx) + ' widthkms 1. linenlam 100'#'lambdarange 5. '#'nlam 100'
print(command)
os.system(command)


plotSpectrum(readSpectrum(), mol=h2o, ilin = lin_idx, fnu = True, dpc = 60.1, jy = True)

#%%
#
# Now make a strongly zoomed-in image
# This time we plot linear intensity, not logarithmic,
# but we saturate the starlight
#
#   radmc3d image lambda 10 incl 60 sizeau 10
#
makeImage( incl = 0, phi = 0, wav = lamb, sizeau = 10, npix = pixel_count)
#
# to get the image
#
imz  = readImage()
dum  = imz.image.shape
nx   = dum[0]
nx2  = nx//2
#vmax = 5e-10
vmax = np.max(imz.image)

fig = plt.figure()
fig.set_size_inches(14,10)
plotImage(imz,au=True,bunit='inu',vmax=vmax,cmap=cm.hot)
plt.show()

#
# Now make sure to have run 
#
#   radmc3d image circ lambda 10 incl 60
#
os.system('radmc3d image circ lambda ' + str(lamb) +' incl 0')
#
# to get the "circular image". Circular images (available only for
# models in spherical coordinates) are pixel arrangements that are
# not as rows and columns, but instead as concentric circles. This 
# can be useful for models with radial coordinates spanning a huge
# range (from very small r to very large r, several orders of 
# magnitude). For circular images we do not need "subpixeling"
# or "zoom ins", because the pixels are arranged in circles, each
# corresponding to a radial coordinate of the model. So this "zoom
# in" is automatic: all scales are represented in one "circular
# image". They are also useful for 1-D models, since 2-D images of
# 1-D spherically symmetric models are overkill: we only need to 
# compute the intensity as a function of impact parameter. So for
# 1-D models the circular images have no subdivision in phi, but
# only intensity as a function of r. That makes it, as a bonus,
# also much faster.
#
# Here we overplot the result of the spherical image with that
# of the (normal) rectangular image.
#
imcirc = readcircimage()
fig = plt.figure()
fig.set_size_inches(14,10)
plt.plot(imcirc.rc,imcirc.image[:,0,0,0],label='spherical image')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('r [cm]')
plt.ylabel(r'$I_\nu [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}\,\mathrm{Hz}^{-1}]$')
plt.plot(im.x[nx2:],im.image[nx2:,nx2,0],'o',label='rectangular image')
#plt.ylim(1e-11, 1e-8)
plt.legend()
plt.show()

#
# You can also show a "heatmap" of the circular image
#
#   The radial span
#   (but note that the r-grid is not fully log, so this is not 100% correct)
#
lgrin  = np.log10(imcirc.rc[1]/au)
lgrout = np.log10(imcirc.rc[-1]/au)
#
#   Make the "heatmap" figure of the 10log of the circular image
#
fig = plt.figure()
fig.set_size_inches(14,10)
plt.imshow(np.log10(imcirc.image[1:,:,0,0]+1e-23),vmin=-16,aspect='auto',cmap=cm.hot,origin='lower',extent=[0,360,lgrin,lgrout])
plt.title(r'$\lambda = 10\,\mu\mathrm{m}$')
plt.xlabel(r'$\phi\; [deg]$')
plt.ylabel(r'$^{10}\log(r)\; [\mathrm{AU}]$')
cbar=plt.colorbar()
cbar.set_label(r'$^{10}\log(I_\nu)\;[\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}\,\mathrm{Hz}^{-1}\,\mathrm{ster}^{-1}\,]$')
plt.show()

