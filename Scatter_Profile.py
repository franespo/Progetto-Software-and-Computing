
import matplotlib.pyplot as plt
import functions as fn
from configparser import ConfigParser


# Choice of wavelength to compute VISibility and normalise particle scattering 
# coefficient
parser = ConfigParser()
parser.read('ScatterProfile_Configuration.ini')

start_wav = parser.getfloat('General_Variables', 'starting_wav',
                         fallback = 0.3)
final_wav = parser.getfloat('General_Variables', 'final_wav',
                         fallback = 1.01)
zena = parser.getfloat('General_Variables', 'zenith_angle',
                         fallback = 40) 
scata = parser.getfloat('General_Variables', 'scattering_angle',
                         fallback = 130) 
maratio = parser.getfloat('General Variables', 'mixing_ratio',
                          fallback = 0)
wav_step = parser.getfloat('General_Variables', 'wavelength_step',
                         fallback = 0.01) 
alpha = parser.getfloat('General_Variables', 'angstrom_coefficient_alpha',
                         fallback = 1) 
g = parser.getfloat('General_Variables', 'Henyey_Greenstein_g_parameter',
                         fallback = 0.85) 
bt = parser.getfloat('General_Variables', 'BT_of_the_sun',
                         fallback = 5800) 
Hv = parser.getfloat('General_Variables', 'Vertical_scale_Height',
                         fallback = 7.7) 
wavelength_norm = parser.getfloat('General_Variables', 'wavelength_norm',
                         fallback = 0.52)
d = parser.getfloat('General_Variables', 'distance_from_black_wall',
                         fallback = 25)

wav, i_wav = fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm)


# FIGURE 1 
tau_vertical = fn.RayOpticaldepth(Hv,wav)
ktot = fn.scatter_coeff(tau_vertical,Hv,wav,i_wav,maratio
                        ,alpha)
k_m = ktot[0]
k_a = ktot[1]   
name_figure = 'Scattering coefficient trend'               
fig = plt.figure()    
plt.plot(wav, k_m, 'b-', label='k-mol.scattering', linewidth = 2)   
plt.plot(wav, k_a, 'r-', label='k-aer.scattering', linewidth = 2)
fig.suptitle(name_figure)
plt.ylabel('Scattering coefficient [1/m]')
plt.xlabel('Wavelength [micron]')
# Definition of legend which is on the plot 
leg = plt.legend()
# Saving figure in the same folder of code
fig.savefig('FiguresNew/fig1.png') 

# FIGURE 2

k_tot = ktot[2]
tau = fn.Opticaldepth(k_tot,Hv,d)
tau_v = tau[0]
tau_h = tau[1]
name_figure = 'Optical Paths'
fig = plt.figure()
plt.plot(wav, tau_v, 'b-', label='Vertical optical path', linewidth = 2)   
plt.plot(wav, tau_h, 'r-', label='Horizontal optical path', linewidth = 2)
fig.suptitle(name_figure)
plt.ylabel('Optical Path')
plt.xlabel('Wavelength [micron]')
# Definition of legend which is on the plot 
leg = plt.legend()
# Saving figure in the same folder of code
fig.savefig('FiguresNew/fig2.png') 

# FIGURE 3

angle = fn.radcos(scata,zena)
transmitt = fn.transmittance(wav,tau_v,tau_h,angle[1])
Ts = transmitt[0]
Th  = transmitt[1]
name_figure = 'Transmittance on vertical and horizontal path (25m)'
fig = plt.figure()
plt.plot(wav, Ts, 'b-', label='Transmittance on vertical path', linewidth = 2)   
plt.plot(wav, Th, 'r-', label='Transimtance on horizontal path', linewidth = 2)
fig.suptitle(name_figure)
plt.ylabel('Transmittance')
plt.xlabel('Wavelength [micron]')
# Definition of legend which is on the plot 
leg = plt.legend()
# Saving figure in the same folder of code
fig.savefig('FiguresNew/fig3.png') 


# FIGURE 4

Stot = fn.ScatteringTot(k_m,k_a,wav,wav_step,g)
ch = Stot[0]
name_figure = 'Scattering diagram'
fig = plt.figure()
# First plot and definition of main properties: semilog y-axis
plt.semilogy(Stot[1], Stot[2], 'b--', label='Ms', linewidth = 2)
# Second plot and definition of main properties: semilog y-axis
plt.semilogy(Stot[1], Stot[3], 'r--', label='As',linewidth = 2)
# Third plot and definition of main properties: semilog y-axis
plt.semilogy(Stot[1], ch[:,0], 'k-', label='Stot(0.30)',linewidth = 3)
# Fourth plot and definition of main properties: semilog y-axis
plt.semilogy(Stot[1], ch[:,1], 'g-', label='Stot(0.64)',linewidth = 3)
# Fifth plot and definition of main properties: semilog y-axis
plt.semilogy(Stot[1], ch[:,2], 'm--', label='Stot(1.0)',linewidth = 3)
fig.suptitle(name_figure)
plt.xlabel('Scattering angle')
plt.ylabel('Scattering diagram')
# Definition of legend which is on the plot 
leg = plt.legend()
# Saving figure in the same folder of code
fig.savefig('FiguresNew/fig4.png')

# FIGURE 5

irrad = fn.irradiance(wav,wav_step,bt,tau_v,angle[1])
name_figure = 'Irradiances'
fig = plt.figure()
plt.plot(wav, irrad[0], 'g-', label='TOA Irradiance', linewidth = 2)
plt.plot(wav, irrad[1], 'b-', label='SL Irradiance',linewidth = 2)
fig.suptitle(name_figure)
plt.xlabel('Wavelength [micron]')
plt.ylabel('Irradiance')
# Definition of legend which is on the plot 
leg = plt.legend()
# Saving figure in the same folder of code
fig.savefig('FiguresNew/fig5.png')

# FIGURE 6
name_figure = 'Normalized Irradiances'
fig = plt.figure()
plt.plot(wav, irrad[2], 'g-', label='TOA Normalized Irradiance', linewidth = 2)
plt.plot(wav, irrad[3], 'b-', label='SL Normalized Irradiance',linewidth = 2)
fig.suptitle(name_figure)
plt.xlabel('Wavelength [micron]')
plt.ylabel('Relative Units')
# Definition of legend which is on the plot 
leg = plt.legend()
# Saving figure in the same folder of code
fig.savefig('FiguresNew/fig6.png')


        

    

    
    
    
    
    
    