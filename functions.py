# Functions used in SingleScatter.py


import math as m
import numpy as np


def wav_matrix(start_wav,wav_step,final_wav,wavelength_norm):
    """ This function computes the wavelength vector
    
        INPUT: start_wav is the starting value of wavelength vector
               wav_step is the step value of wavelength vector
               final_wav is the final value of wavelength
               
        OUTPUT: wavelength vector       
    """
    
    # Wavelength vector
    wav = np.arange(start_wav,final_wav,wav_step)
    
    # Conversion of wavelength vector in wavelength matrix                                      
    wav = wav.T
    
    i_wav = np.round((wavelength_norm-start_wav)/wav_step)
        
    return wav, i_wav

def radcos(scata,zena):
    """ This function computes the cosines of scattering and sun zenith angles
    
        INPUT: scata is the scattering angle expressed in degrees
               zena is the sun scattering angle expressed in degrees
               
        OUTPUT: cosines of scattering and sun zenith angles expressed in 
                radians
    """
    
    # Cosine of scattering angle expressed in radians 
    ang_mu = np.cos(scata*np.pi/180)            
    # Cosine of sun zenith angle expressed in radians 
    ang_mus = np.cos(zena*np.pi/180)
    
    return ang_mu, ang_mus


def RayOpticaldepth(Hv,wav):
    """ This function computes the Rayleigh optical depth matrix
        
        INPUT: Hv is the Vertical scale height
               wav is the wavelength vector
               
        OUTPUT: tau_vertical is the Rayleigh optical depth
    """
       
    # Definition of Rayleigh optical depth matrix  
    tau_vertical = np.zeros((np.size(wav)))
    
    for i in range(np.size(wav)):
        # Definition of coeff to compute tau_vertical
        coeff = 3.916 + 0.074 * wav[i] + 0.05 / wav[i]
        # Rayleigh optical depth
        tau_vertical[i] = 8.38e-3 / wav[i]**coeff
        
    return tau_vertical

def scatter_coeff(tau_vertical,Hv,wav,i_wav,maratio, alpha):
    """ This function computes the scattering coefficients
        
        INPUT: tau_vertical is the Rayleigh optical depth
               Hv is the Vertical scale height
            
        OUTPUT: k_mol_scattering, k_aer_scattering and k_tot_scattering
                are respectively the molecular, aerosol and total scattering
                coefficients
    """
    # Rayleigh molecular scattering coefficient 
    k_mol_scattering = tau_vertical / Hv                          
    # Rayleigh aerosol scattering coefficient
    k_aer_scattering = tau_vertical / Hv  
    # Total scattering coefficient                             
    k_tot_scattering = k_mol_scattering + k_aer_scattering
    # Normalize aerosol scattering coefficient using maratio at specific 
    # wavelength 
        
    k_aernorm = k_mol_scattering[int(i_wav)]*maratio/wav[int(i_wav)
                ]**(-alpha)
    for i in range(np.size(wav)):
        k_aer_scattering[i] = k_aernorm * wav[i]**(-alpha)
    
    return k_mol_scattering, k_aer_scattering, k_tot_scattering

def Opticaldepth(k_tot_scattering,Hv):
    """ This function computes the vertical and horizontal optical depth (for
        a path of length d)
    
        INPUT: k_tot_scattering is the total scattering coefficient
               Hv is the Vertical scale height
               d is the distance from observer to black wall expressed in m
               
        OUTPUT: tau_v and tau_h are respectively the vertical and horizontal 
                optical depth
    """
    # Total vertical optical depth
    tau_v = Hv*k_tot_scattering      
    # Distance from observer to black wall [m]   
    d = 25            
    # Total horizontal optical depth for a path of length d
    tau_h = d*k_tot_scattering     
        
    return tau_v, tau_h

def dir_scatteringTot(tau_vertical,Hv,ang_mu,g):
    """ This function computes the total scattering direction
    
    Parameters: 
        dir_mol is scattering direction of molecules (dir_mol)
        dir_aer is scattering direction of aerosols (dir_aer)
        
        Return:
            Total scattering direction expressed by a single value """
    # Rayleigh molecular scattering coefficient
    k_mol_scattering = tau_vertical / Hv
    # Rayleigh aerosol scattering coefficient                 
    k_aer_scattering = k_mol_scattering
    # Scattering direction of molecules
    dir_scattering_mol = 0.750*(1 + ang_mu*ang_mu) 
    # Scattering direction of aerosols D_aer is assumed independent of wa           
    dir_scattering_aerosol = (1-g**2)/(1+g**2-2*g*ang_mu)**(3/2)
    # Total Scattering direction                        
    dir_tot = (dir_scattering_mol*k_mol_scattering + 
             dir_scattering_aerosol*k_aer_scattering
             )/(k_mol_scattering+k_aer_scattering)
    
    return dir_tot[0], dir_tot

# Definition of parameters for Planck function
# Planck function


def emittance(wav,bt):
    """ This function computes the monochromatic emittance of a black body at 
        a specific temperature (brightness temperature) and depends on 
        wavelength.
    
    Parameters:
        wav is the wavelength express as a vector
        bt is the brightness temperature
        i is the index of wavelength vector
        
        Constant:
            c0 is the speed of light in vaccum, 2.997e+8 
            h is Planck constanr, 6.625e-34
            k_ is Boltzmann constant, 1.38e-23 
            nn is the refravtive index of the medium, 1
            
            Raise:
                brightness temperature cannot be negative
                
                Return:
                intensity is the monochromatic emittance """  
                
    c0 = 2.997e+8                 # m/s speed of light in vaccum
    h = 6.625e-34                 # J.s Planck constant 
    kb = 1.38e-23                 # T/K Boltzmann constant
    nn = 1        
    
    if bt <= 0:
        raise ValueError("bt must to be greater than zero")
    intensity = (((2*np.pi*h*(c0**2))/(wav**5))*(
                        1/(np.exp((h*c0)/(kb*wav*bt))-1)))*nn**(-2)
    return intensity


def transmittance(wav,tau_v,tau_h,ang_mus):
    """ definizione funzione
    """
    
    # Definition of Transmittance matrix on vertical path
    Ts = np.zeros((np.size(wav)))
    # Definition of Transmittance on horizontal path matrix
    Th = np.zeros((np.size(wav)))
       
    for i  in range((np.size(wav))):
        # Transmittance on vertical path
        Ts[i] = np.exp(-tau_v[i]/ang_mus)
        # Transmittance on horizontal path
        Th[i] = np.exp(-tau_h[i])
    return Ts, Th

def irradiance(wav,wav_step,bt,tau_v,ang_mus):
    """ descrizione funzione
    """
    # Definition of Irradiance matrix reaching Top of Atmosphere (TOA)
    E0 = np.zeros((np.size(wav)))
    # Definition of Transmittance matrix on vertical path
    Ts = np.zeros((np.size(wav)))
    # Definition of Irradiance matrix reaching Surface Layer (SL)
    E0S = np.zeros((np.size(wav)))
    
    # Coefficient to compute the irradiance on a horizontal path
    B = 1.0e-06       
                 
    for i  in range(np.size(wav)):
        # Irradiance reaching TOA
        E0[i] = emittance(wav[i],bt) * B * ang_mus
        # Transmittance on vertical path
        Ts[i] = np.exp(-tau_v[i]/ang_mus)
        # Irradiance reaching SL
        E0S[i] = E0[i]*Ts[i]
    # Normalization of E0 and E0S with E0N and E0SN
    E0max = np.max(E0)
    E0Smax = np.max(E0S)
    E0N = E0/E0max
    E0SN = E0S/E0Smax
    return E0, E0S, E0N, E0SN

def ScatteringTot(k_mol_scattering,k_aer_scattering,wav,wav_step,g):
    """ This function computes the Total Scattering which depends on Molecular
    Scattering and Aerosols Scattering.
    
    Parameters:
        Mol_scattering is the molecular scattering defined by a matrix which 
        depends on cosine of angle expressed in radians
        Aer_scattering is the aerosols scattering defined by a matrix which 
        depends on cosine of angle expressed in radians and on 
        Henyey-Greenstein g parameter
                
                Return:
                Stot is total scattering which depends on Rayileigh molecular 
                scattering coefficient and Rayileigh aerosol scattering 
                coefficient """               
    # Definition of angles range (0-180) for x-axis with different steps
    # step for range angle 0-10
    sxstep = 0.5
    # Definiton of angle vector 1
    ang1 = np.arange(0,10.1,sxstep)
    # step for range angle 10-180
    dxstep = 5
    # Definition of angle vector 2
    ang2 = np.arange(11,182,dxstep)
    # Union of angle vector 1 and angle vector 2 which represents angle range
    ang = np.concatenate((ang1,ang2),axis=0)
    # Definition of raw size matrix angle
    na = np.size(ang)
    # Definition of angle matrix formed by the cosine of the values ​​of the angle 
    # vector ang expressed in radians
    ang_mua = np.zeros((na,1))
    for i in range(na):
        # New angle matrix
        ang_mua[i] = m.cos(ang[i]*(m.pi/180))
    
    # Definition of a list
    was = [1, m.floor(np.size(wav)/2), np.size(wav)]
    # For scattering direction plot
    # Definition of molecular scattering matrix
    Ms = np.zeros((na,1))
    # Definition of aerosols scattering matrix
    As = np.zeros((na,1))
    # Definition of column size total scattering matrix with dimension of was
    sa = np.size(was)
    # Definition of Total Scattering matrix
    Stot = np.zeros((na,sa))

    for j in range(na):
        if j < 0:
            raise ValueError("index must be greater than zero")
        # Molecular scattering
        Ms[j] = 0.750*(1 + ang_mua[j]*ang_mua[j]) 
        # Aerosol scattering     
        As[j] =(1-g**2) / (1+g**2-2*g*ang_mua[j])**(3/2)   
        for k in range(sa):
            if k < 0:
                raise ValueError("index must be greater than zero")
            # Total Scattering
            Stot[j,k] = (Ms[j]*k_mol_scattering[k]+ 
                        As[j]*k_aer_scattering[k])/(k_mol_scattering[k]+
                          k_aer_scattering[k])
    return Stot, ang, Ms, As
