# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 11:15:55 2020

@author: Francesco
"""
import SingleScatter as sin
import math as m

# Definition of total scattering direction

def dir_scatteringTot(dir_mol,dir_aer):
    # Rayleigh molecular scattering coefficient
    k_mol_scattering = sin.tauf_vertical / sin.Hv
    # Rayleigh aerosol scattering coefficient                 
    k_aer_scattering = k_mol_scattering
    # Total Scattering direction                        
    D_tot = (D_mol*k_mol_scattering + 
             D_aer*k_aer_scattering)/(k_mol_scattering+k_aer_scattering)  
    return D_tot[0]

# Scattering direction of molecules
D_mol = 0.750*(1 + sin.ang_mu*sin.ang_mu) 
# Scattering direction of aerosols D_aer is assumed independent of wa           
D_aer = (1-sin.g**2)/(1+sin.g**2-2*sin.g*sin.ang_mu)**(3/2)
# Total Scattering direction  
D_tot = dir_scatteringTot(D_mol,D_aer)       

# Definition of parameters for Planck function
# Planck function
c0 = 2.997e+8                 # m/s speed of light in vaccum
h = 6.625e-34                 # J.s Planck constant 
k_ = 1.38e-23                 # T/K Boltzmann constant
nn = 1                        # Refravtive index of the medium

def planck(wav,bt,i):
    if bt <= 0:
        raise ValueError(
                "bt (brightness temperature) must to be greater than zero")
    a = ((2*m.pi*h*(c0**2))/(wav**5))
    b = (1/(m.exp((h*c0)/(k_*wav*(i+1)*bt))-1))
    intensity = (a*b)*nn**(-2)
    return intensity


def ScatteringTot(Mol_scattering,Aer_scattering,j,k):             
        for j in range(sin.na):
            if j < 0:
                raise ValueError("index must be greater than zero")
            # Molecular scattering
            sin.Ms[j] = 0.750*(1 + sin.ang_mua[j]*sin.ang_mua[j]) 
            # Aerosol scattering     
            sin.As[j] =(1-sin.g**2) / (1+sin.g**2
                  -2*sin.g*sin.ang_mua[j])**(3/2)   
            for k in range(sin.sa):
                if k < 0:
                    raise ValueError("index must be greater than zero")
                # Total Scattering
                sin.Stot[j,k] = (sin.Ms[j]*sin.k_mol_scattering[k]
                + sin.As[j]*sin.k_aer_scattering[k])/(sin.k_mol_scattering[k]
                +sin.k_aer_scattering[k])
        return sin.Stot
