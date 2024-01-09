# This file is part of Hydrasim, a small piece of software intended to 
# simulate observations of TW Hydra.

# Hydrasim is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# Hydrasim is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with Hydrasim. If not, see <https://www.gnu.org/licenses/>. 

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 10:15:25 2023

@author: cassie
"""

import casatasks
import casatools
import os
import numpy as np
import matplotlib.pyplot as plt
#import simalma
#import numpy as np

twhya_coord = 'J2000 11h01m51.9054s -34d42m17.0316s'

file = "oh2o1600.fits"

project_name = "oh2o_sim_hires27"

vp = casatools.vpmanager()



# casatasks.imhead(file)

# ia = casatools.image()
# ia.open(infile=file)


# axesLength = ia.shape()
# # Divide the first two elements of axesLength by 2.
# center_pixel = [ x / 2.0 for x in axesLength[:2] ]
# # Feed center_pixel to ia.toworld and and save the RA and Dec to ra_radians and dec_radians
# (ra_radians, dec_radians) = ia.toworld( center_pixel )['numeric'][:2]

# qa = casatools.quanta()
# ra_hms  = qa.formxxx(str(ra_radians)+"rad",format='hms',prec=5)
# dec_dms = qa.formxxx(str(dec_radians)+"rad",format='dms',prec=5)

# mycs = ia.coordsys()

# ia.regrid(outfile = "co_regrid2.fits", 
#           shape = [1000,1000, 40, 1],
#           axes = [0,1],
#           csys = mycs.torecord(),
#           overwrite = True)

# ia.close()

#%%

basecfg = "scifi_proto.cfg"

dir = '/home/cassie/casa/casa-6.5.2-26-py3.8/data/alma/simmos/'

xyscale = 12
xoffset = xyscale * 2000/50
yoffset = xyscale * 28000/50
zout = 0
diamout = 4

outconfig = "scifi.cfg"
cfg = "scifi"

f_in = open(dir + basecfg, 'r')
f_out = open(outconfig, 'w')

for line in f_in:
    if (line[0] == ('#')):
        f_out.write(line)
        print(line)
    else:
        fields = line.split()
        x = fields[0]
        y = fields[1]
        z = fields[2]
        diam = fields[3]
        pad = fields[4]
        
        x = str(float(x) * xyscale + xoffset)
        y = str(float(y) * xyscale + yoffset)
        z = str(zout)
        diam = str(diamout)
        
        outdata = [x,y,z,diam,pad,"\n"]
        f_out.write(" ".join(outdata))
    
f_out.close()
f_in.close()

#vp.setpbairy("SCIFI", dishdiam = diamout, blockagediam = 0)
vp.setpbairy("SCIFI", dishdiam = "2m",blockagediam = "0m")

# hpbw 
# vp.setpbairy("SCIFI", reffreq = "557GHz", 

scivp = vp.getvp("SCIFI",freq=500e9)

#%% do some integration time estimation

kb = 1.381e-23 # boltzmann content

Jy_convert = 1e-26 # conversion from W/m/m/Hz to Jy

noiselev = 1e-3 #janskies per beam

bw = 3.71e6 # bandwidth in Hz
Tsys = 25 # system temperature in kelvin
efficiency = 0.8 * 0.9 * 0.85 # efficiencies multiplied together
n_ant = 10 # number of antennas

A = np.pi * (diamout/2)**2 # area of antenna in m^2

SEFD = Tsys*(2*kb/(efficiency*A)) # SEFD in SI units

SEFD_Jy = SEFD/Jy_convert # convert to Jy/beam

tint = (SEFD_Jy/noiselev)**2 / (n_ant*(n_ant-1)*bw) # integration time required in seconds

tint_h = tint / 3600 # integration time in hours

tint_runs = tint_h/14 # number of runs required

print(str(tint_h))



#%%

# casatasks.simobserve(
#     project = project_name,
#     skymodel = file,
#     setpointings       =  True,
#     direction          =  twhya_coord,
#     indirection        =  twhya_coord,
#     #incell = "0.005arcsec",
#     #mapsize            =  "0.76arcsec",
#     obsmode            =  "int",
#     totaltime          =  "14h",
#     antennalist        =  cfg + ".cfg",
#     thermalnoise       =  '')


# apply_noise = True
# if (apply_noise):
#     sm = casatools.simulator()
#     sm.openfromms(project_name + "/" + project_name + "." + cfg + ".ms")
#     sm.setseed(seed=int(11215))
#     sm.setnoise(
#         mode = 'tsys-manual',
#         trx = float(275),
#         tau = float(0.0), 
#         rxtype = int(2))  #rxtype 1 is 2SB, 2 is DSB
#     sm.corrupt()
    
#     sm.done()
    
#%%

numnights = tint_runs
night_vises = []
vis_prefix = project_name + "/" +project_name + "."


#%%
for night in range(0, 1):
    print(night)
    
    
    cfg_temp= cfg
    
    if not os.path.exists(cfg_temp+".cfg"):
        os.link(outconfig, cfg_temp+".cfg")
    
    night_vises.append(vis_prefix + cfg_temp + ".ms")
    
    casatasks.simobserve(
            project = project_name,
            skymodel = file,
            setpointings       =  True,
            direction          =  twhya_coord,
            indirection        =  twhya_coord,
            #incell = "0.005arcsec",
            #mapsize            =  "0.76arcsec",
            obsmode            =  "int",
            totaltime          =  "14h",
            antennalist        =  cfg + ".cfg",
            thermalnoise       =  '')
    
    
    apply_noise = True
    if (apply_noise):
            sm = casatools.simulator()
            sm.openfromms(project_name + "/" + project_name + "." + cfg + ".ms")
            sm.setseed(seed=int(10000 + night))
            sm.setnoise(
                mode = 'tsys-manual',
                # here we reduce the receiver temperature by a factor equal to the 
                # sqrt of the number of observing passes conducted to emulate the noise
                # reduction from performing that many observing passes
                
                # this is a bit hacky but seems like the easiest way to do it.
                trx = float(Tsys)/np.sqrt(float(numnights)),
                tcmb= float(2.73)/np.sqrt(float(numnights)),
                tground=float(2.73)/np.sqrt(float(numnights)),
                tatmos=float(12),
                tau = float(0.0), 
                rxtype = int(2))  #rxtype 1 is 2SB, 2 is DSB
            sm.corrupt()
            
            sm.done()
    #%%
    

# night_vises = []
# vis_prefix = project_name + "/" +project_name + "."

# for night in range(0, numnights):
#     cfg_temp="scifi" + str(night);
#     night_vises.append(vis_prefix + cfg_temp + ".ms")
    
modelimage = vis_prefix + cfg + ".skymodel.flat" 
  
#%%  
# note: to run tclean from scratch, delete scifi.model
# otherwise, tclean will simply apply more clean iterations to existing model
# casatasks.tclean(
#         vis = project_name + "/" +project_name + "." + cfg + ".ms",
#         imagename= project_name + "/" +project_name + "." + cfg,
#         #startmodel=modelimage,
#         imsize = [600,600],
#         cell="0.005arcsec",
#         niter = 0,
#         threshold = "1e-4Jy",
#         weighting = "natural"
#         )

cleanprior = "cleanprior"

ia = casatools.image()
ia.open(modelimage) 

im2 = ia.rebin(outfile = cleanprior, bin=[160,160],overwrite=True)
im2.done()
ia.close()
#%%

imsize_pixels = 1000

casatasks.simanalyze(
    project = project_name,
    image = True,
    #modelimage=cleanprior,
    vis = "$project." + cfg + ".ms",# + "," + project_name + "." + cfg + ".ms",
    imsize = imsize_pixels,
    niter = 20000,
    threshold = "1e-7Jy",
    weighting = "natural",
    analyze = True,
    showuv = False,
    showresidual = True,  
    showconvolved = True,
    showclean = True,
    graphics = "both",
    verbose = True,
    overwrite = True)

#%% 
spacing = 10
means = []
radii = np.arange (0.01,imsize_pixels/2, spacing)

for i in radii:
    region = "annulus[[" + str(imsize_pixels/2) + "pix, " + str(imsize_pixels/2) + "pix], [" + str(i) + "pix, " + str(i + spacing) + "pix]]"
    
    print(region)
    
    statout = casatasks.imstat(
        imagename= project_name + "/" + project_name + ".scifi.image.flat", 
        region=region,
        )
    
    means.append(statout['mean'][0])
    
headout = casatasks.imhead(
    imagename= project_name + "/" + project_name + ".scifi.image.flat")

radtodeg = 180/np.pi
degtoarcsec = 3600

pixsize = abs(headout['incr'][0]) * radtodeg * degtoarcsec

plt.plot(radii*pixsize, means)
plt.title("Radially Averaged Intensity")
plt.xlabel("Distance from disk center (arcsec)")
plt.ylabel("Average Intensity (Jy/beam)")
plt.xlim(0,max(radii*pixsize))
plt.ylim(-0.01, 0.03)
plt.show()
    

#%%

print(str(Tsys) + "K")
print(str(diamout) + "m") 
print(str(tint_h/24) + " days")
print(str())

# #%%
# casatasks.simalma(
#     project = project_name,
#     overwrite=True,
#     skymodel = file,
#     #indirection="J2000 23h59m59.96s -34d59m59.50s",
#     #incell="0.1arcsec",
#     #inbright="0.004",
#     #incenter="330.076GHz",
#     #inwidth="50MHz",
#     antennalist=["alma.cycle8.3.cfg"],
#     totaltime="1800s",
#     #tpnant=2,
#     #tptime="7200s",
#     pwv=0.6,
#     #mapsize="1arcmin",
#     niter=100,
#     dryrun=False,
#     imsize=[1024,1024])