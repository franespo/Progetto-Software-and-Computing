# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 10:43:26 2020

@author: Francesco
"""

# This code is implemented to investigate the single scattering solutions and
# the "colours" of the light reaching the observer from an horizontal path. The
# atmosphere may contain molecules and aerosol. 

import numpy as np
import math

# Definition of parameters that don't change in the code

g = 0.85    # Henyey-Greenstein g parameter
            
bt = 5800   # Brightness Temperature (BT) of the Sun

# Different cases are valuated expressed by job='A#', where # is the number
# of case analyzed

# Definition of parameters that change in the code

job = ('A0')
zena = 40                # The Sun zenith angle
scata = 130              # The scattering angle
maratio = 0              # The ratio of molecules to aerosol scattering at 0.56
alpha = 1                # The Angstrom coefficient


# Choice of wavelength to compute VISibility and normalise particle scattering 
# coefficient

wd = 0.01                # Step in wavelength vector
iwa = 27                 # wa(iwa)= 0.52 microns

wa = np.arange(0.3,1.0,wd)          # Wavelength vector from VIS to NIR
n = np.size(wa)                     # Number of steps             
wa = np.arange(0.3,1.0,wd).reshape((n,1)) # Wavelength matrix
mu = math.cos(scata*math.pi/180)    # Conversion from degrees to radians
mus = math.cos(zena*math.pi/180)    # Conversion from degrees to radians

# How to compute Rayleigh optical depth tauf[i] for a vertical path through 
# the atmosphere with airmass=1

tauf = np.zeros((n,1))

for i in range(n):
    barg = 3.916 + 0.074 * wa[i] + 0.05 / wa[i]
    tauf[i] = 8.38e-3 / wa[i]**barg
    
# Vertical scale height, Hv, is the same for molecular scattering coefficient 
# ksm and aerosol scattering coefficient ksa useful to compute them

Hv = 7.7 
ksm = tauf / Hv                       # Rayleigh molecular scattering coefficient    
ksa = tauf / Hv                       # Rayleigh aerosol scattering coefficient        

# Normalize aerosol scattering coefficient using maratio

k_aernorm=ksm[iwa]*maratio/wa[iwa]**(-alpha)






       
      
        

    
    

    



