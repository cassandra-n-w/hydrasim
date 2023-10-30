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
#import simalma
#import numpy as np

twhya_coord = 'J2000 11h01m51.9054s -34d42m17.0316s'

file = "oh2o1600.fits"

project_name = "oh2o_sim_hires13"

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

xyscale = 10
xoffset = xyscale * 2000/50
yoffset = xyscale * 28000/50
zout = 0
diamout = 2

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



#%%

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
    sm.setseed(seed=int(11215))
    sm.setnoise(
        mode = 'tsys-manual',
        trx = float(275),
        tau = float(0.0), 
        rxtype = int(2))  #rxtype 1 is 2SB, 2 is DSB
    sm.corrupt()
    
    sm.done()
    
#%%
cfg="scifi2"

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
        sm.setseed(seed=int(11216))
        sm.setnoise(
            mode = 'tsys-manual',
            trx = float(275),
            tau = float(0.0), 
            rxtype = int(2))  #rxtype 1 is 2SB, 2 is DSB
        sm.corrupt()
        
        sm.done()
#%%

casatasks.tclean(
        vis = [project_name + "/" +project_name + "." + "scifi" + ".ms", project_name + "/" +project_name + "." + cfg + ".ms"],
        imagename= project_name + "/" +project_name + "." + cfg,
        imsize = 600,
        cell="0.005arcsec",
        niter = 0000,
        threshold = "1e-7Jy",
        weighting = "natural"
        )

casatasks.simanalyze(
    project = project_name,
    image = False,
    #modelimage = file,
    vis = "$project." + cfg + ".ms",# + "," + project_name + "." + cfg + ".ms",
    imsize = [600, 600],
    niter = 10000,
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