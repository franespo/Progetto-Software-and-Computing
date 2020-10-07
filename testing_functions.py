#!/usr/bin/python3

# Testing functions.py

import numpy as np
import functions as fn
import pytest
from hypothesis import strategies as st
from hypothesis import settings
from hypothesis import given

# Test for wav_matrix
 
@given(st.floats(0.1,4),st.floats(0.5,1),st.floats(10,15),st.floats(5,9.99))
@settings(max_examples = 25)

def test_wav_matrix(start_wav,wav_step,final_wav,wavelength_norm):

    wav, i_wav = fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm)
    
    #check if output is positive
    assert(i_wav >= 0)
    assert(len(wav[wav < 0]) ==0)
    
    #check if i_wav values is a correct and valid index
    assert(wav[i_wav]>0)
    #check if i_wav is conteined in wav vector
    assert(i_wav <= (len(wav)-1))
    
    #check if a ValueError arise if start_wav is smaller than 0
    with pytest.raises(ValueError):
       fn.wav_matrix(-1,0.01,5,2)
        
    #check if ValueError arise if final_wav is smaller than start_wav
    with pytest.raises(ValueError):
       fn.wav_matrix(0.5,0.1,0.4,0.3)
       

# Test for radcos

@given(st.floats(0,180),st.floats(0,90))
@settings(max_examples = 25)

def test_radcos(scata,zena):
	
    ang_mu, ang_mus = fn.radcos(scata,zena)
    
    #check if cosine value is correct
    assert(ang_mu >= -1 and ang_mu <= 1) 
    assert(ang_mus >= -1 and ang_mus <= 1)

    #check if a ValueError arise if scata is smaller than 0
    with pytest.raises(ValueError):
       fn.radcos(-5,10)
        
    #check if a ValueError arise if zena is smaller than 0
    with pytest.raises(ValueError):
       fn.radcos(10,-5)
      
    #check if a ValueError arise if scata is greater than 90
    with pytest.raises(ValueError):
       fn.radcos(200,10)
        
    #check if a ValueError arise if zena is bigger than 180
    with pytest.raises(ValueError):
       fn.radcos(10,200)

# Test for RayOpticaldepth
       
@given(st.floats(0.1,4),st.floats(0.5,1),st.floats(10,15),st.floats(5,9.99))
@settings(max_examples = 20)
def test_RayOpticaldepth(start_wav,wav_step,final_wav,wavelength_norm):

    wav, _ = fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm)
    
    
    tau_vertical = fn.RayOpticaldepth(wav)
    #check if the output length is the correct
    assert(len(tau_vertical) == len(wav))
        
        
# Test for scatter_coeff
       
@given(st.floats(0.1,100),st.floats(0.1,4),st.floats(0.5,1),st.floats(10,15),
	st.floats(5,9.99),st.floats(0,50),st.floats(0,20))
@settings(max_examples = 20)

def test_scatter_coeff(Hv,start_wav,wav_step,final_wav,wavelength_norm,maratio, alpha):

    wav, i_wav = fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm)
    tau_vertical = fn.RayOpticaldepth(wav)
    k_mol_scattering, k_aer_scattering, k_tot_scattering = fn.scatter_coeff(tau_vertical,
    							Hv,wav,i_wav,maratio, alpha)
    
    #check if a ValueError arise if Hv is smaller or equal to 0
    with pytest.raises(ValueError):    
    	 fn.scatter_coeff(tau_vertical,-1,wav,i_wav,0, 1)
    	 7
    #check if the output length is correct
    assert(len(k_mol_scattering) == len(tau_vertical))
    

# Test for Opticaldepth
    
@given(st.floats(0.1,100),st.floats(0,100),st.floats(0.1,4),st.floats(0.5,1),
	st.floats(10,15),st.floats(5,9.99))
@settings(max_examples = 20)

def test_Opticaldepth(Hv,d,start_wav,wav_step,final_wav,wavelength_norm):

    wav, _ = fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm)
    k_tot_scattering = np.random.rand(len(wav))
    tau_h, tau_v = fn.Opticaldepth(k_tot_scattering,Hv,d)
    
    #check if the output length is correct
    assert(len(tau_v) == len(wav))
    
    #check if a ValueError arise if d is 0
    with pytest.raises(ValueError):
        fn.Opticaldepth(k_tot_scattering,Hv,-5)    

# Test for dir_scatteringTot

@given(st.floats(0.1,4),st.floats(0.5,1),st.floats(10,15),st.floats(5,9.99),
	st.floats(0.1,100),st.floats(-1,1))
@settings(max_examples = 20)
        
def test_dir_scatteringTot(start_wav,wav_step,final_wav,wavelength_norm,Hv,ang_mu):

    wav, _ = fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm)
    tau_vertical = fn.RayOpticaldepth(wav)
    
    dir_tot = fn.dir_scatteringTot(tau_vertical,Hv,ang_mu)
       
    #check if output do not contain negative elements
    assert(dir_tot[0] > 0)
   
# Test for emittance

@given(st.floats(0.1,4),st.floats(0.5,1),st.floats(10,15),st.floats(5,9.99),
	st.floats(1,8000))
@settings(max_examples = 20)

def test_emittance(start_wav,wav_step,final_wav,wavelength_norm,bt):
   
   wav, _ = fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm)
   intensity = fn.emittance(wav,bt)
   
   #check if output is positive
   assert(len(intensity[intensity <= 0]) == 0)
   
   #testing with a negative value, the value error is detected
   with pytest.raises(ValueError):
       fn.emittance(wav,-1)
       
   #check if the output length is correct
   assert(len(intensity) == len(wav))
   
# Test for transmittance

@given(st.floats(0.1,4),st.floats(0.5,1),st.floats(10,15),st.floats(5,9.99),
	st.floats(0,200),st.floats(0,200),st.floats(0.1,1))
@settings(max_examples = 20)

def test_transmittance(start_wav,wav_step,final_wav,wavelength_norm,tau_v,tau_h,ang_mus):
   
   wav, _ = fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm) 
   tau_h = np.random.rand(len(wav))
   tau_v = np.random.rand(len(wav))
   
   Ts, Th = fn.transmittance(wav,tau_v,tau_h,ang_mus)
   
   
   #check if output is positive
   assert(len(Ts[Ts <=0]) == 0)
   assert(len(Th[Th <=0]) == 0)
   
   #check if length of th is equal to length of wav
   assert(len(Ts) == len(wav))
   
   #check if length of ts is equal to length of wav
   assert(len(Th) == len(wav))
  

# Test for irradiance
   
@given(st.floats(0.1,4),st.floats(0.5,1),st.floats(10,15),st.floats(5,9.99),
	st.floats(1,8000),st.floats(0,100), st.floats(0.1,1))
@settings(max_examples = 20)

def test_irradiance(start_wav,wav_step,final_wav,wavelength_norm,bt,tau_v,ang_mus):
   
   wav, _ = fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm)
   tau_v = np.random.rand(len(wav))
   E0, E0S, E0N, E0SN = fn.irradiance(wav,wav_step,bt,tau_v,ang_mus)
   
   #check if output is positive
   assert(len(E0[E0 <=0]) == 0)
   assert(len(E0S[E0S <=0]) == 0)
   assert(len(E0N[E0N <=0]) == 0)
   assert(len(E0SN[E0SN <=0]) == 0)
    
   #check if length of E0 is equal to length of wav
   assert(len(E0) == len(wav))
   #check if length of E0 is equal to length of wav
   assert(len(E0S) == len(wav))
   #check if length of E0 is equal to length of wav
   assert(len(E0N) == len(wav))
   #check if length of E0 is equal to length of wav
   assert(len(E0SN) == len(wav))
   
   
# Test for ScatteringTot
   
@given(st.floats(0,100),st.floats(0,100),st.floats(0.1,4),st.floats(0.5,1),
	st.floats(10,15),st.floats(5,9.99))
@settings(max_examples = 20)

def test_ScatteringTot(k_mol_scattering,k_aer_scattering,start_wav,wav_step,
			final_wav,wavelength_norm):
    
    wav, _ = fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm)
    k_mol_scattering = np.random.rand(len(wav))
    k_aer_scattering = np.random.rand(len(wav))
    
    Stot, ang, Ms, As = fn.ScatteringTot(k_mol_scattering,k_aer_scattering,wav,wav_step)
    
    #check if output is positive
    assert(len(Stot[Stot <=0]) == 0)
    assert(len(ang[ang < 0]) == 0)
    assert(len(Ms[Ms <=0]) == 0)
    assert(len(As[As <=0]) == 0)

   

if __name__ == '__main__':
    pass
   
   
