# This file is part of Hydrasim, a small piece of software intended to 
# simulate observations of TW Hydra.

# Hydrasim is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# Hydrasim is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with Hydrasim. If not, see <https://www.gnu.org/licenses/>. 


#%% Import libraries and define functions
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import scipy 
from scipy import interpolate
import math
import matplotlib.pyplot as plt


"""
Created on Wed Mar 23 10:51:18 2022

@author: cassie
"""

def load_data_as_dic(filepath, comments='!', returnOriginalKeys=False):
    """load iter_000*.dat into a dictionary 
       from Fujun Du's pylib"""

    data = np.loadtxt(filepath, comments=comments)

    ftmp = open(filepath, 'r')
    str_comment = ftmp.readline()[1:].split()
    ftmp.close()

    dic = {}
    for i in range(len(str_comment)):
        dic.update({str_comment[i]: data[:, i]})

    del data

    if returnOriginalKeys:
        return str_comment, dic
    else:
        return dic
    
def resample(dic, index, rhoc, thetac):
    
    rhrhc, ttc = np.meshgrid(rhoc, thetac, indexing = 'ij');
    
    rc = np.cos(ttc)*rhrhc
    zc = np.sin(ttc)*rhrhc
    
    rcinterp = 0.5 * (dic['rmin'] + dic['rmax'])
    zcinterp = 0.5 * (dic['zmin'] + dic['zmax'])
    datainterp = np.array(dic[index]);
    
    interpolant = scipy.interpolate.interp2d(rcinterp, zcinterp, datainterp);
    rc = rc.ravel(order='F')   
    zc = zc.ravel(order='F')
    
    resampled = interpolant(rc, zc);
    
    return resampled[0,:]

def resample2(dic, index, rhoc, thetac):
    
    #rhrhc, ttc = np.meshgrid(rhoc, thetac, indexing = 'ij');
    
    #rc = np.cos(ttc)*rhrhc
    #zc = np.sin(ttc)*rhrhc
    
    rcinterp = 0.5 * (dic['rmin'] + dic['rmax'])
    zcinterp = 0.5 * (dic['zmin'] + dic['zmax'])
    thetacinterp = np.arctan2(rcinterp, zcinterp)
    rhocinterp = np.sqrt(rcinterp**2.0 +zcinterp**2.0)
    
    datainterp = np.array(dic[index]);
    fill = np.min(np.log(datainterp)) - 1
    print(fill)
    interp2 = (np.pi/2 - thetacinterp)
    interpolant = scipy.interpolate.interp2d(np.log(rhocinterp), interp2, np.log(datainterp), fill_value=fill);
    #rc = rc.ravel(order='F')   
    #zc = zc.ravel(order='F')
    
    resampled = np.exp(interpolant(np.log(rhoc), (np.pi/2 - thetac)));
    
    resampled = np.clip(resampled, 0, np.max(datainterp))
    
    return resampled
    

# THIS IS THE GOOD ONE, USE THIS ONE
def resample3(dic, index, rhoc, thetac, fillmod = 0):
    
    rhrhc, ttc = np.meshgrid(rhoc, thetac, indexing = 'ij');
    
    #rc = np.cos(ttc)*rhrhc
    #zc = np.sin(ttc)*rhrhc
    
    rcinterp = 0.5 * (dic['rmin'] + dic['rmax'])
    zcinterp = 0.5 * (dic['zmin'] + dic['zmax'])
    thetacinterp = np.arctan2(rcinterp, zcinterp)
    rhocinterp = np.sqrt(rcinterp**2.0 +zcinterp**2.0)
    
    #datainterp = np.array(dic[index]);
    datainterp = index;
    fill = np.min(np.log(np.abs(datainterp))) - fillmod
    print(fill)
    zcinterp_mod = (np.pi/2 - thetacinterp)
    #interpolant = scipy.interpolate.interp2d(np.log(rhocinterp), interp2, np.log(datainterp), fill_value=fill);
    
    points_table = (np.log(rhocinterp), zcinterp_mod)
    #points_sample = (np.log(rhoc), (np.pi/2 - thetac))
    points_sample = (np.log(rhrhc), (np.pi/2 - ttc))
    
    resampled = scipy.interpolate.griddata(points_table, np.log(np.abs(datainterp)), points_sample, fill_value=fill)
    
    
    
    #resampled = np.clip(resampled, 0, np.max(datainterp))
    
    # do some tricky magic to make samples inside the disk but outside the sampling boundary
    # take on the value of whatever is above them
    
    # flip the array so that the center comes first, and find the boolean value
    # of whether or not the object is inbounds
    oob = np.flip(resampled == fill, axis=1)
    
    
    # perform a sequential and operation using cumprod to find all INNER DISK LOCATIONS ONLY
    # that are out of bounds
    oob = np.cumprod(oob, axis=1)
    
    # find the first inbound element
    first_inbound = np.minimum(np.sum(oob,axis=1), np.size(resampled, axis=1) - 1)
    
    
    first_inbound_val = np.flip(resampled, axis=1)
    first_inbound_val = first_inbound_val[np.arange(np.size(resampled, axis=0)), first_inbound]
    
    #first_inbound_val = np.tile(first_inbound_val, (1,np.size(resampled, axis=1)))
    first_inbound_val = np.reshape(first_inbound_val,(first_inbound_val.size, 1))
    resampled = np.where(oob, first_inbound_val, np.flip(resampled, axis=1))
    
    resampled = np.flip(resampled, axis=1)
    
    
    resampled = (np.exp(resampled));
    
    return resampled

def resample4(dic, index, rhoc, thetac):
    
    rhrhc, ttc = np.meshgrid(rhoc, thetac, indexing = 'ij');
    
    #rc = np.cos(ttc)*rhrhc
    #zc = np.sin(ttc)*rhrhc
    
    rcinterp = 0.5 * (dic['rmin'] + dic['rmax'])
    zcinterp = 0.5 * (dic['zmin'] + dic['zmax'])
    thetacinterp = np.arctan2(rcinterp, zcinterp)
    rhocinterp = np.sqrt(rcinterp**2.0 +zcinterp**2.0)
    
    datainterp = np.array(dic[index]);
    fill = np.min(np.log(datainterp)) - 1
    print(fill)
    zcinterp_mod = (np.pi/2 - thetacinterp)
    #interpolant = scipy.interpolate.interp2d(np.log(rhocinterp), interp2, np.log(datainterp), fill_value=fill);
    
    points_table = (np.log(rhocinterp), zcinterp_mod)
    #points_sample = (np.log(rhoc), (np.pi/2 - thetac))
    points_sample = (np.log(rhrhc), (np.pi/2 - ttc))
    
    resampled = scipy.interpolate.griddata(points_table, np.log(datainterp), points_sample, fill_value=fill)
    
    resampled = (np.exp(resampled));
    
    #resampled = np.clip(resampled, 0, np.max(datainterp))
    
    return resampled


#%% Load data from file
    

path = "/home/cassie/hydrasim/model_data/iter_0001.dat"

dic = load_data_as_dic(path)

#%% preprocess coordinate data to convert cylindrical to spherical

zmin = dic['zmin']
zmax = dic['zmax']
rmin = dic['rmin']
rmax = dic['rmax']
rcinterp = 0.5 * (dic['rmin'] + dic['rmax'])
zcinterp = 0.5 * (dic['zmin'] + dic['zmax'])


thetamin = np.arctan2(rmin,zmin)
thetamax = np.arctan2(rmax,zmax)
rhomin = np.sqrt(rmin**2.0 + zmin**2.0)
rhomax = np.sqrt(rmax**2.0 + zmax**2.0)





#%% Setup problem parameters and resample data
# following code based on run_spher2d example
# Some natural constants
#
au  = 1.49598e13     # Astronomical Unit       [cm]
pc  = 3.08572e18     # Parsec                  [cm]
ms  = 1.98892e33     # Solar mass              [g]
ts  = 5.78e3         # Solar temperature       [K]
ls  = 3.8525e33      # Solar luminosity        [erg/s]
rs  = 6.96e10        # Solar radius            [cm]
G = 6.6743e-08        # gravitational constant  [dyn cm^2 / g^2]
#
# Monte Carlo parameters
#
nphot    = 100000
#
# Grid parameters
#
nx       = 100
ny       = 1500        # One-sided only. So the "real" value is twice this.
nz       = 1

# sample theta more closely to the disk plane
#thetai = np.flip((1 - np.logspace(-2, 0,ny)) * np.pi/2);
thetai = np.linspace(np.pi/2. - np.pi/3.8, np.pi/2, ny+1)
thetai[-1] = np.pi/2.;
thetai[0] = 0

# sample rho logarithmically also
rhoi = np.logspace(np.log10(np.min(rhomax)), np.log10(np.max(rhomin) * 1.5), nx+1)




#
# Model parameters
# Irrelevant,comes from example model
# rin      = 5*au
# rout     = 100*au
# zmaxr    = 0.7e0
# rho0     = 1e-16 * 10000
# prho     = -2.e0
# hpr      = 0.1e0
#
# Star parameters
# I took these from wikipedia
# star values for TW Hya
mstar    = 0.8*ms
rstar    = 1.11*rs
tstar    = 4000 #Kelvin
pstar    = [0.,0.,0.]

#rstar = rstar * 0.0005

#
# Make the coordinates
#
xi       = rhoi
yi       = thetai
zi       = np.array([0.,math.pi*2])
xc       = 0.5e0 * ( xi[:-1] + xi[1:] )
yc       = 0.5e0 * ( yi[:-1] + yi[1:] )



#
# Make the dust density model
#
# rr,tt    = np.meshgrid(xc,yc,indexing='ij')
# zzr      = math.pi/2.0 - tt
# rhod2     = rho0 * (rr/au)**prho
# rhod2     = rhod2 * np.exp(-0.50*(zzr/hpr)**2)

md_cell = dic['md_cell'];

T_gas_dic = dic['Tgas'];

T_dust_dic = dic['Tdust']

mg_cell = dic['mg_cell']

#CO_ratio = dic['mg_cell']
CO_ratio = dic['CO']

H2O_ratio = dic['H2O']

n_gas_dic = dic['n_gas']

mass_CO = 28.01 #grams/mol
mass_H = 1.00784
mass_HD = 3.02204
mass_H2 = 2 * mass_H
avogadro = 6.022140e23 #molecules / mole

molecular_mass_CO = mass_CO/avogadro;



#calculate volume of each cell in cm^3
vol_cell = np.pi * (zmax - zmin) * (rmax**2 - rmin**2) * au**3

dust_mass_density_zr = md_cell/vol_cell
gas_density = mg_cell/vol_cell

# calculate rho_CO from total gas mass per cell
# note! average mass of gas per molecule is 1.4 * hydrogen mass, because of
# high quantities of helium

# THIS IS WRONG DO NOT FUCKING USE IT, it calculates the mass density (correctly)
# but we want the NUMBER DENSITY for molecular lines
rho_CO_unsampled1 = CO_ratio * (mass_CO / (mass_H * 1.4)) * gas_density;

# calculate rho_CO from gas number density

rho_CO_unsampled2 = CO_ratio * n_gas_dic;
rho_H2O_unsampled = H2O_ratio * n_gas_dic;

#1e-10 * 
rho_co1 = resample3(dic, rho_CO_unsampled1, xc, yc, 60)
rho_co2 = resample3(dic, rho_CO_unsampled2, xc, yc, 60)

# multiply by 0.5, assuming 50:50 ortho/para water ratio
rho_oh2o = 0.5 * resample3(dic, rho_H2O_unsampled, xc, yc, 60)

test_ratio = rho_CO_unsampled1/rho_CO_unsampled2
rho_co = rho_co2

rhod = resample3(dic, dust_mass_density_zr, xc, yc, 30)
T_gas = resample3(dic, T_gas_dic, xc, yc)
n_gas = resample3(dic, n_gas_dic, xc, yc, 30)
T_dust = resample3(dic, T_dust_dic, xc, yc)

mass_tot = resample3(dic, md_cell + mg_cell, xc, yc, 10)

mass_tot_r = np.sum(mass_tot, axis=1)
mass_tot_within_r = np.cumsum(mass_tot_r) + mstar

veloc = np.sqrt(G * mass_tot_within_r / (xc * au))
veloc = np.transpose(np.tile(veloc, (ny, 1)))

dust = dust_mass_density_zr



#
# Write the wavelength_micron.inp file
#
lam1     = 0.1e0
lam2     = 7.0e0
lam3     = 25.0e0
lam4     = 1.0e4
n12      = 20
n23      = 100
n34      = 30
lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
lam      = np.concatenate([lam12,lam23,lam34])
nlam     = lam.size

#%% perform some debug plotting

#plt.style.use('_mpl-gallery-nogrid')

xcgrid, ycgrid = np.meshgrid(xc, yc);

rc = np.sin(ycgrid)*xcgrid
zc = np.cos(ycgrid)*xcgrid

#rcm, zcm = np.meshgrid(rc, zc)

interp = (np.transpose(T_gas)).flatten()
un_interp = T_gas_dic

#interp = (np.transpose(n_gas)).flatten()
#un_interp = n_gas_dic

rcm = rc.flatten()
zcm = zc.flatten()

colormap = "jet";

#levels = np.linspace(np.log10(np.min(interp)), np.log10(np.max(interp)), 30);
levels = np.linspace(np.log10(6), np.log10(2000), 200)

fig, (ax1, ax2) = plt.subplots(1,2)
fig.set_size_inches(14,6)
#ax1.plot(np.log10(rcinterp), np.log10(zcinterp), 'b.')
#ax1.tricontour(np.log10(rcinterp), np.log10(zcinterp), np.log10(un_interp), levels=levels)
#ax1.plot((rcinterp), (zcinterp), 'b.')
tc = ax1.tricontourf((rcinterp), (zcinterp), np.log10(un_interp), levels=levels, cmap = colormap)
ax1.set(xlim = (0, 140), ylim = (0, 140))

ax1.set_ylabel("Z coordinate (AU)")
ax1.set_xlabel("R coordinate (AU)")
ax1.set_title("Temperature (log10 K) of Protoplanetary Disk (Du Model)")



lrcm = np.log10(rcm)
lzcm = np.log10(zcm)
lrhodp = np.maximum(np.log10(interp), np.array(-70))
#ax2.plot(lrcm, lzcm, 'b.')
#ax2.tricontour(np.log10(rcm), np.log10(zcm), lrhodp, levels=levels)    
#ax2.plot(10**lrcm, 10**lzcm, 'b.')
tc = ax2.tricontourf((rcm), (zcm), lrhodp, levels=levels, cmap = colormap)  
fig.colorbar(tc)

ax2.set_ylabel("Z coordinate (AU)")
ax2.set_xlabel("R coordinate (AU)")
ax2.set_title("Temperature (log10 K) of Protoplanetary Disk (Interpolated)")

ax2.set(xlim=ax1.get_xlim(),ylim = ax1.get_ylim())
    
plt.show()

#%% do more plotting

un_interp = T_gas_dic

#levels = np.linspace(np.log10(np.min(interp)), np.log10(np.max(interp)), 30);
levels = np.linspace(np.log10(6), np.log10(2000), 200)

fig, ax1 = plt.subplots(1,1)
fig.set_size_inches(8,6)
#ax1.plot(np.log10(rcinterp), np.log10(zcinterp), 'b.')
#ax1.tricontour(np.log10(rcinterp), np.log10(zcinterp), np.log10(un_interp), levels=levels)
#ax1.plot((rcinterp), (zcinterp), 'b.')
tc = ax1.tricontourf((rcinterp), (zcinterp), np.log10(un_interp), levels=levels, cmap = colormap)
ax1.set(xlim = (0, 190), ylim = (0, 170))

ax1.set_ylabel("Z coordinate (AU)")
ax1.set_xlabel("R coordinate (AU)")
ax1.set_title("Gas Temperature (log10 K) (Du Model)")

fig.colorbar(tc)

    
plt.show()

#%% do more plotting

un_interp = T_dust_dic

#levels = np.linspace(np.log10(np.min(interp)), np.log10(np.max(interp)), 30);
levels = np.linspace(np.log10(6), np.log10(200), 200)

fig, ax1 = plt.subplots(1,1)
fig.set_size_inches(8,6)
#ax1.plot(np.log10(rcinterp), np.log10(zcinterp), 'b.')
#ax1.tricontour(np.log10(rcinterp), np.log10(zcinterp), np.log10(un_interp), levels=levels)
#ax1.plot((rcinterp), (zcinterp), 'b.')
tc = ax1.tricontourf((rcinterp), (zcinterp), np.log10(un_interp), levels=levels, cmap = colormap)
ax1.set(xlim = (0, 190), ylim = (0, 170))

ax1.set_ylabel("Z coordinate (AU)")
ax1.set_xlabel("R coordinate (AU)")
ax1.set_title("Dust Temperature (log10 K) (Du Model)")

fig.colorbar(tc)

    
plt.show()

#%% do more plotting

un_interp = abs(rho_H2O_unsampled)

levels = np.linspace(np.log10(np.min(un_interp)), np.log10(np.max(un_interp)), 30);
#levels = np.linspace(np.log10(6), np.log10(2000), 200)


fig, ax1 = plt.subplots(1,1)
fig.set_size_inches(8,6)
#ax1.plot(np.log10(rcinterp), np.log10(zcinterp), 'b.')
#ax1.tricontour(np.log10(rcinterp), np.log10(zcinterp), np.log10(un_interp), levels=levels)
#ax1.plot((rcinterp), (zcinterp), 'b.')
tc = ax1.tricontourf((rcinterp), (zcinterp), np.log10(un_interp), levels=levels, cmap = colormap)
ax1.set(xlim = (0, 190), ylim = (0, 170))

ax1.set_ylabel("Z coordinate (AU)")
ax1.set_xlabel("R coordinate (AU)")
ax1.set_title("H2O Density (log10 1/cm^3) (Du Model)")

fig.colorbar(tc)

    
plt.show()

#%% do more plotting

un_interp = abs(rho_CO_unsampled2)

levels = np.linspace(np.log10(np.min(un_interp)), np.log10(np.max(un_interp)), 30);
#levels = np.linspace(np.log10(6), np.log10(2000), 200)


fig, ax1 = plt.subplots(1,1)
fig.set_size_inches(8,6)
#ax1.plot(np.log10(rcinterp), np.log10(zcinterp), 'b.')
#ax1.tricontour(np.log10(rcinterp), np.log10(zcinterp), np.log10(un_interp), levels=levels)
#ax1.plot((rcinterp), (zcinterp), 'b.')
tc = ax1.tricontourf((rcinterp), (zcinterp), np.log10(un_interp), levels=levels, cmap = colormap)
ax1.set(xlim = (0, 190), ylim = (0, 170))

ax1.set_ylabel("Z coordinate (AU)")
ax1.set_xlabel("R coordinate (AU)")
ax1.set_title("CO Density (log10 1/cm^3) (Du Model)")

fig.colorbar(tc)

    
plt.show()


#%% Write Output Files
#
# Write the wavelength file
#
xi_au = xi * au
rho_out = rhod


# write the lines file
with open('lines.inp', 'w+') as f:
    f.write('2\n') # mode, set to 2
    f.write('1\n') # number of molecular/atomic species
    #f.write('co leiden 0 0 0')
    f.write('oh2o leiden 0 0 0')
    
with open('wavelength_micron.inp','w+') as f:
    f.write('%d\n'%(nlam))
    for value in lam:
        f.write('%13.6e\n'%(value))
        
        
#
#
# Write the stars.inp file
#
with open('stars.inp','w+') as f:
    f.write('2\n')
    f.write('1 %d\n\n'%(nlam))
    f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar,mstar,pstar[0],pstar[1],pstar[2]))
    for value in lam:
        f.write('%13.6e\n'%(value))
    f.write('\n%13.6e\n'%(-tstar))
#
# Write the grid file
#
with open('amr_grid.inp','w+') as f:
    f.write('1\n')                       # iformat
    f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
    f.write('100\n')                     # Coordinate system
    f.write('0\n')                       # gridinfo
    f.write('1 1 0\n')                   # Include x,y,z coordinate
    f.write('%d %d %d\n'%(nx,ny,nz))     # Size of grid
    for value in xi_au:
        f.write('%13.6e\n'%(value))      # X coordinates (cell walls)
    for value in yi:
        f.write('%13.8e\n'%(value))      # Y coordinates (cell walls)
    for value in zi:
        f.write('%13.6e\n'%(value))      # Z coordinates (cell walls)
#
# Write the density file
#
with open('dust_density.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nx*ny*nz))           # Nr of cells
    f.write('1\n')                       # Nr of dust species
    data = rho_out.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')

#write the CO density file
with open('numberdens_co.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nx*ny*nz))           # Nr of cells
    data = rho_co.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')
    
#write the H2O density file
with open('numberdens_oh2o.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nx*ny*nz))           # Nr of cells
    data = rho_oh2o.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')
    
#write the H2O density file
with open('numberdens_ph2o.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nx*ny*nz))           # Nr of cells
    data = rho_oh2o.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')        
    
#write the gas temperature file
with open('gas_temperature.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nx*ny*nz))           # Nr of cells
    data = T_gas.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')
    
#write the dust temperature file
with open('dust_temperature.dat','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nx*ny*nz))           # Nr of cells
    f.write('1\n')                       # Nr of dust species
    data = T_dust.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')
    
#write the dust temperature file
with open('gas_velocity.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nx*ny*nz))           # Nr of cells
    data = veloc.ravel(order='F')        # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="0.000000e+00  0.000000e+00 %13.6e")
    f.write('\n')


#
# Dust opacity control file
#
with open('dustopac.inp','w+') as f:
    f.write('2               Format number of this file\n')
    f.write('1               Nr of dust species\n')
    f.write('============================================================================\n')
    f.write('1               Way in which this dust species is read\n')
    f.write('0               0=Thermal grain\n')
    f.write('silicate        Extension of name of dustkappa_***.inp file\n')
    f.write('----------------------------------------------------------------------------\n')
#
# Write the radmc3d.inp control file
#
with open('radmc3d.inp','w+') as f:
    f.write('nphot = %d\n'%(nphot))
    f.write('scattering_mode_max = 0\n')   # Put this to 1 for isotropic scattering
    f.write('temp1 = 1.5e5\n') #increasing max temperature for like reasons
    f.write('lines_partition_temp1 = 1.5e5\n') #still doing the increase max temperature
