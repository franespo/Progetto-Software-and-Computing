# This code is implemented to investigate the single scattering solutions and
# the "colours" of the light reaching the observer from an horizontal path. The
# atmosphere may contain molecules and aerosol. 

import numpy as np
import math as m
import matplotlib.pyplot as plt

# Definition of parameters that don't change in the code
# Henyey-Greenstein g parameter
g = 0.85    

# Brightness Temperature (BT) of the Sun           
bt = 5800   

# Different cases are valuated
# Some suggestions on which values ​​to use within the code to obtain 
# different cases:

# zena=40, scata=130, maratio=0, alpha=1, i_wav = 26
# zena=40, scata=130, maratio=1, alpha=1, i_wav = 26
# zena=40, scata=130, maratio=1, alpha=3, i_wav = 26
# zena=40, scata=130, maratio=10, alpha=1, i_wav = 26
# zena=70, scata=160, maratio=10, alpha=1, i_wav = 26
# zena=40, scata= 50, maratio=10, alpha=1, i_wav = 26
    

# Definition of parameters that change in the code

while(True):
    try:
        # The Sun zenith angle
        zena = int(input('Insert Sun zenith angle "zena": ')) 
        if(zena >=0 and zena <=90):
            break
        else:
            print("zena must have a value from 0 to 90. Try again")
    except ValueError:
        print("You must insert numbers not words/letters")
                    
while(True):
    try:
        # The scattering angle
        scata = int(input('Insert scattering angle "scata": ')) 
        if(scata >=0 and scata <=180):
            break
        else:
            print("scata must have a value from 0 to 180. Try again")
    except ValueError:
        print("You must insert numbers not words/letters")
        
while(True):
    try:
        # The ratio of molecules to aerosol scattering at 0.56 micron(i_wav=26)
        maratio = int(input(
                'Insert ratio of molecules to aerosol scattering "maratio": ')) 
        if(maratio >=0 and maratio <=10):
            break
        else:
            print("maratio must have a value from 0 to 10. Try again")
    except ValueError:
        print("You must insert numbers not words/letters")

while(True):
    try:
        # The Angstrom coefficient
        alpha = int(input('Insert Angstrom coefficient "alpha": '))
        if(alpha >=0):
            break
        else:
            print("alpha must have a value > 0. Try again")
    except ValueError:
        print("You must insert numbers not words/letters")        

# Choice of wavelength to compute VISibility and normalise particle scattering 
# coefficient

# Step in wavelength vector
wav_step = 0.01         
# Wavelength vector from VIS (0.3) to NIR (1.01)
wav = np.arange(0.3,1.01,wav_step)
# Number of steps n in wavelength vector wav
n = np.size(wav)
# Conversion of wavelength vector in wavelength matrix                                      
wav = np.arange(0.3,1.01,wav_step).reshape((n,1))
# Cosine of scattering angle expressed in radians 
ang_mu = m.cos(scata*m.pi/180)            
# Cosine of sun zenith angle expressed in radians 
ang_mus = m.cos(zena*m.pi/180)

# How to compute Rayleigh optical depth tauf[i] for a vertical path through 
# the atmosphere with airmass = 1

# Definition of Rayleigh optical depth matrix
tauf_vertical = np.zeros((n,1))
for i in range(n):
    # Definition of coeff to compute tauf_vertical
    coeff = 3.916 + 0.074 * wav[i] + 0.05 / wav[i]
    # Rayleigh optical depth
    tauf_vertical[i] = 8.38e-3 / wav[i]**coeff
    
# Vertical scale height, Hv, is the same for molecular scattering coefficient 
# k_mol_scattering and aerosol scattering coefficient k_aer_scattering

# Vertical Scale height
Hv = 7.7                                                
# Rayleigh molecular scattering coefficient 
k_mol_scattering = tauf_vertical / Hv                          
# Rayleigh aerosol scattering coefficient
k_aer_scattering = tauf_vertical / Hv                               

# This input is essential to choose which wavelength you want to put through
# the i-th position of the wavelength matrix 
while(True):
    try:
        # i_wav is the i-th position of wavelength matrix: i_wav = 0 
        # corrispondsto 0.3 micron, i_wav = 1 corrisponds to 0.31 micron, 
        # i_wav = 2 corrisponds to 0.32 micron and so on and so forth
        i_wav =int(input("Insert the value of i_wav that is the i-position of wavelength vector and must have a value between 0 and n. For the case at 0.56 micron, it is recommanded 26: ")) 
        if(i_wav >=0 and i_wav <=n):
            break
        else:
            print("i_wav must have a value between 0 and n. Try again")
    except ValueError:
        print("You must insert numbers not words/letters")
        
 # Normalize aerosol scattering coefficient using maratio at specific 
 # wavelength       
k_aernorm = k_mol_scattering[i_wav]*maratio/wav[i_wav]**(-alpha)
for i in range(n):
    k_aer_scattering[i] = k_aernorm * wav[i]**(-alpha)

# How to compute optical depth for a vertical path and for a horizontal path

# Distance from observer to black wall [m]   
d = 25
# Total scattering coefficient                             
k_tot_scattering = k_mol_scattering + k_aer_scattering 
# Total vertical optical depth for a path of length d
tau_v = Hv*k_tot_scattering        
# Total horizontal optical depth for a path of length d
tau_h = d*k_tot_scattering         

# Definition of total scattering direction
    
def dir_scatteringTot(dir_mol,dir_aer):
    """ This function computes the total scattering direction.
    
    Parameters: 
        dir_mol is scattering direction of molecules (dir_mol)
        dir_aer is scattering direction of aerosols (dir_aer)
        
        Return:
            Total scattering direction expressed by a single value """
            
    # Rayleigh molecular scattering coefficient       
    k_mol_scattering = tauf_vertical / Hv
    # Rayleigh aerosol scattering coefficient                        
    k_aer_scattering = k_mol_scattering
    # Total Scattering direction  [1/sr]                     
    D_tot = (D_mol*k_mol_scattering + 
             D_aer*k_aer_scattering)/(k_mol_scattering+k_aer_scattering)  
    return D_tot[0]

# Scattering direction of molecules
D_mol = 0.750*(1 + ang_mu*ang_mu)
# Scattering direction of aerosols D_aer is assumed independent of wa           
D_aer = (1-g**2)/(1+g**2-2*g*ang_mu)**(3/2)
# Total Scattering direction 
D_tot = dir_scatteringTot(D_mol,D_aer)       
# Connection with VISibility
vis = (3.912 / (k_mol_scattering[i_wav]+k_aer_scattering[i_wav]))

# Definition of parameters for Planck function
# Planck function

# speed of light in vacuum
c0 = 2.997e+8
# Planck constant                 
h = 6.625e-34
# Boltzmann constant 
k_ = 1.38e-23
# Refravtive index of medium                  
nn = 1                        

def planck(wav,bt,i):
    """ This function computes the monochromatic emittance of a black body at a
    specific temperature (brightness temperature) and depends on wavelength.
    
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

    if bt <= 0:
        raise ValueError("bt must to be greater than zero")
    a = ((2*m.pi*h*(c0**2))/(wav**5))
    b = (1/(m.exp((h*c0)/(k_*wav*(i+1)*bt))-1))
    intensity = (a*b)*nn**(-2)
    return intensity

# Definition of Transmittance matrix on vertical path
Ts = np.zeros((n,1))
# Definition of Irradiance matrix reaching Top of Atmosphere (TOA)          
E0 = np.zeros((n,1))
# Definition of Irradiance matrix reaching Surface Layer (SL)        
E0S = np.zeros((n,1))
# Definition of Monochromatic emittance matrix        
intensity = np.zeros((n,1))
# Coefficient to compute the irradiance on a horizontal path
B = 1.0e-06                     
   
for i  in range(n):
    intensity[i] = planck(wav[i],bt,i)
    # Transmittance on vertical path
    Ts[i] = m.exp(-tau_v[i]/ang_mus)
    # Irradiance reaching TOA   
    E0[i] = intensity[i] * B * ang_mus
    # Irradiance reaching SL    
    E0S[i] = E0[i]*Ts[i]
    
# Definition of Transmittance on horizontal path matrix
Th = np.zeros((n,1))
# Definition of Radiance matrix
L = np.zeros((n,1))
                     
for i in range(n):
    # Transmittance on horizontal path
    Th[i] = m.exp(-tau_h[i])
    # Radiance            
    L[i] = E0S[i] * D_tot * (1-Th[i])
    
# Normalization of L, E0 and E0S with LN, E0N and E0SN
LN = np.zeros((n,1))
Lmax = np.max(L)
E0max = np.max(E0)
E0Smax = np.max(E0S)
LN = L/Lmax
E0N = E0/E0max
E0SN = E0S/E0Smax

# Definition of titles on the graphics using already defined values 
# The function str(A) converts a numeric array into a string array 
# that represents the numbers

st0 = str(zena) 
st1 = str(scata)
s = " "
ss = ", "
st00 = s.join(['sun zenith angle =', st0])
st11 = s.join(['scattering angle =', st1])
str0 = ss.join([st00,st11])

st0 = str(d)
st1 = str(maratio)
st2 = str(vis)
st3 = str(alpha)
st00 = s.join(['d =',st0])
st11 = s.join(['M/A =',st1])
st22 = s.join(['vis =',st2])
st33 = s.join(['alpha =',st3])
str1 = ss.join([st00,st11,st22,st33])

# Figures
# Creation of figure 1 with dimension 20x20
fig = plt.figure(figsize=(20,20))
# Creation of figure with subplots, in this case the final figure will have
# two plots: definition of first subplot
ax = fig.add_subplot(2,1,1)
# First plot and definition of main properties
ax.plot(wav, k_mol_scattering, 'b-', label='k_mol_scattering', linewidth = 3)
# Second plot and definition of main properties
ax.plot(wav, k_aer_scattering, 'r-', label='k_aer_scattering',linewidth = 3)
# Definition of subplot title
plt.title(str0)
# Definition of y-axis name
plt.ylabel('Scattering coeff. [1/m]')
# Definition of legend which is on the plot 
leg = ax.legend()

# Definition of second subplot
ax = fig.add_subplot(2,1,2)
# First plot and definition of main properties
ax.plot(wav, tau_v, 'b-', label='tau_v', linewidth = 3)
# Second plot and definition of main properties
ax.plot(wav, tau_h, 'r-', label='tau_h',linewidth = 3)
# Definition of subplot title
plt.title(str1)
# Definition of x-axis name
plt.xlabel('Wavelength [micron]')
# Definition of y-axis name
plt.ylabel('Optical Path')
# Definition of legend which is on the plot
leg = ax.legend()
# Saving figure in the same folder of code
fig.savefig('Figures/fig1.png')

# Creation of figure 2 with dimension 20x20
fig = plt.figure(figsize=(20,20))
# Creation of figure with subplots, in this case the final figure will have
# two plots: definition of first subplot
ax = fig.add_subplot(2,1,1)
# First plot and definition of main properties
ax.plot(wav, Ts, 'b-', label='Transmittance on vertical path', linewidth = 3)
# Second plot and definition of main properties
ax.plot(wav, Th, 'r-', label='Trasmittance on horizontal path',linewidth = 3)
# Definition of subplot title
plt.title(str1)
# Definition of x-axis name
plt.xlabel('Wavelength [micron]')
# Definition of y-axis name
plt.ylabel('Transmittance')
# Definition of legend which is on the plot 
leg = ax.legend()

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
was = [1, m.floor(n/2), n]
# For scattering direction plot
# Definition of molecular scattering matrix
Ms = np.zeros((na,1))
# Definition of aerosols scattering matrix
As = np.zeros((na,1))
# Definition of column size total scattering matrix with dimension of was
sa = np.size(was)
# Definition of Total Scattering matrix
Stot = np.zeros((na,sa))

def ScatteringTot(Mol_scattering,Aer_scattering):
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
    
    for j in range(na):
        # Molecular scattering
        Ms[j] = 0.750*(1 + ang_mua[j]*ang_mua[j]) 
        # Aerosol scattering     
        As[j] =(1-g**2) / (1+g**2-2*g*ang_mua[j])**(3/2)   
        for k in range(sa):
            # Total Scattering
            Stot[j,k] = (Ms[j]*k_mol_scattering[k]
        + As[j]*k_aer_scattering[k])/(k_mol_scattering[k]+k_aer_scattering[k])
    return Stot

Stot = ScatteringTot(Ms,As)

# Definition of second subplot
ax = fig.add_subplot(2,1,2)
# First plot and definition of main properties: semilog y-axis
ax.semilogy(ang, Ms, 'b--', label='Ms', linewidth = 3)
# Second plot and definition of main properties: semilog y-axis
ax.semilogy(ang, As, 'r--', label='As',linewidth = 3)
# Third plot and definition of main properties: semilog y-axis
ax.semilogy(ang, Stot[:,0], 'k-', label='Stot(0.30)',linewidth = 3)
# Fourth plot and definition of main properties: semilog y-axis
ax.semilogy(ang, Stot[:,1], 'g-', label='Stot(0.64)',linewidth = 3)
# Fifth plot and definition of main properties: semilog y-axis
ax.semilogy(ang, Stot[:,2], 'm--', label='Stot(1.0)',linewidth = 3)
# Definition of x-axis name
plt.xlabel('Scattering angle')
# Definition of y-axis name
plt.ylabel('Scattering diagram')
# Definition of legend which is on the plot 
leg = ax.legend()
# Saving figure in the same folder of code
fig.savefig('Figures/fig2.png')

# Creation of figure 3 with dimension 20x20
fig = plt.figure(figsize=(20,20))
# Creation of figure with subplots, in this case the final figure will have
# two plots: definition of first subplot
ax = fig.add_subplot(2,1,1)
# First plot and definition of main properties
ax.plot(wav, E0, 'g-', label='E0', linewidth = 3)
# Second plot and definition of main properties
ax.plot(wav, E0S, 'b-', label='E0S',linewidth = 3)
# Third plot and definition of main properties
ax.plot(wav, L, 'r-', label='L',linewidth = 3)
# Definition of subplot title
plt.title(str1)
# Definition of x-axis name
plt.xlabel('Wavelength [micron]')
# Definition of y-axis name
plt.ylabel('Irradiance')
# Definition of legend which is on the plot
leg = ax.legend()

# Definition of second subplot
ax = fig.add_subplot(2,1,2)
# First plot and definition of main properties
ax.plot(wav, E0N, 'g-', label='E0N', linewidth = 3)
# Second plot and definition of main properties
ax.plot(wav, E0SN, 'b-', label='E0SN',linewidth = 3)
# Third plot and definition of main properties
ax.plot(wav, LN, 'r-', label='LN',linewidth = 3)
# Definition of x-axis name
plt.xlabel('Wavelength [micron]')
# Definition of y-axis name
plt.ylabel('Relative units')
# Definition of legend which is on the plot
leg = ax.legend()
# Saving figure in the same folder of code
fig.savefig('Figures/fig3.png')


