#!/usr/bin/python3

# Testing functions.py

import numpy as np
import functions as fn
import pytest
from hypothesis import strategies as st
from hypothesis import settings
from hypothesis import given

# Test for function 'wav_matrix'
 
@given(st.floats(0.1,4),st.floats(0.5,1),st.floats(10,15),st.floats(5,9.99))
@settings(max_examples = 10)

def test_wav_matrix(start_wav,wav_step,final_wav,wavelength_norm):

    wav, i_wav = fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm)
    
    #check if output of function is positive
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
       

# Test for function 'radcos'

@given(st.floats(0,180),st.floats(0,90))
@settings(max_examples = 10)

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

# Test for function 'RayOpticaldepth'
       
@given(st.floats(0.1,4),st.floats(0.5,1),st.floats(10,15),st.floats(5,9.99))
@settings(max_examples = 5)

def test_RayOpticaldepth(start_wav,wav_step,final_wav,wavelength_norm):

    wav, _ = fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm)
    tau_vertical = fn.RayOpticaldepth(wav)
    
    #check if the output length is the correct
    assert(len(tau_vertical) == len(wav))
    
    #check if the output is positive
    assert(len(tau_vertical[tau_vertical <=0]) == 0) 
    
 
# Test for function 'scatter_coeff'
       
@given(st.floats(0.1,100),st.floats(0.1,4),st.floats(0.5,1),st.floats(10,15),
       st.floats(5,9.99),st.floats(1,50))
@settings(max_examples = 20)

def test_scatter_coeff(Hv,start_wav,wav_step,final_wav,wavelength_norm,
                       maratio):

    wav, i_wav = fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm)
    tau_vertical = fn.RayOpticaldepth(wav)
    k_mol_scattering, k_aer_scattering, k_tot_scattering = fn.scatter_coeff(
            tau_vertical,Hv,wav,i_wav,maratio)
    
    #check if a ValueError arise if Hv is smaller or equal to 0
    with pytest.raises(ValueError):    
    	 fn.scatter_coeff(tau_vertical,-5,wav,i_wav,0)
    #check if a ValueError arise if maratio is negative
    with pytest.raises(ValueError):
    	fn.scatter_coeff(tau_vertical,Hv,wav,i_wav,-10)

    #check if maratio is zero so k_aer_scattering have all zero elements
    if maratio == 0:
       assert(k_aer_scattering[k_aer_scattering == 0])
       	 
    #check if the output length is correct
    assert(len(k_mol_scattering) == len(tau_vertical))
    #check if the output length is correct
    assert(len(k_aer_scattering) == len(tau_vertical))
    #check if the output length is correct
    assert(len(k_tot_scattering) == len(tau_vertical))

# Test for function 'Opticaldepth'
    
@given(st.floats(0.1,100),st.floats(1,100),st.floats(0.1,4),st.floats(0.5,1),
       st.floats(10,15),st.floats(5,9.99))
@settings(max_examples = 20)

def test_Opticaldepth(Hv,d,start_wav,wav_step,final_wav,wavelength_norm):

    wav, _ = fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm)
    k_tot_scattering = np.random.rand(len(wav))
    tau_h, tau_v = fn.Opticaldepth(k_tot_scattering,Hv,d)
    
    #check if the output length is correct
    assert(len(tau_v) == len(k_tot_scattering))
    
    #check if k_tot_scattering is composed by zero elements
    #so tau_v and tau_h also have all zero elements
    assert(len(tau_h[tau_h ==0]) == len(k_tot_scattering[k_tot_scattering ==0])) 
    assert(len(tau_v[tau_v ==0]) == len(k_tot_scattering[k_tot_scattering ==0])) 
    
    #check if a ValueError arise if d is zero
    with pytest.raises(ValueError):
        fn.Opticaldepth(k_tot_scattering,Hv,0)    
    #check if a ValueError arise if Hv is smaller or equal to 0
    with pytest.raises(ValueError):    
    	 fn.Opticaldepth(k_tot_scattering,-10,d)
    	 
# Test for function 'dir_scatteringTot'

@given(st.floats(0.1,4),st.floats(0.5,1),st.floats(10,15),st.floats(5,9.99),
       st.floats(0.1,100),st.floats(-1,1))
@settings(max_examples = 10)
        
def test_dir_scatteringTot(start_wav,wav_step,final_wav,wavelength_norm,
                           Hv,ang_mu):

    wav, _ = fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm)
    tau_vertical = fn.RayOpticaldepth(wav)
    dir_tot = fn.dir_scatteringTot(tau_vertical,Hv,ang_mu)
              
    #check if the output do not contain negative elements
    assert(len(dir_tot[dir_tot < 0]) == 0)   
    #check if the output has the same length of input tau_vertical
    assert(len(dir_tot) == len(tau_vertical))
    
    #check if cosine has value out of the range [-1,1]
    with pytest.raises(ValueError):
         fn.dir_scatteringTot(tau_vertical,Hv,-5)
            
# Test for function 'emittance'

@given(st.floats(0.1,4),st.floats(0.5,1),st.floats(10,15),st.floats(5,9.99),
       st.floats(1,8000))
@settings(max_examples = 10)

def test_emittance(start_wav,wav_step,final_wav,wavelength_norm,bt):
   
   wav, _ = fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm)
   intensity = fn.emittance(wav,bt)
   
   #check if the output is positive
   assert(len(intensity[intensity <= 0]) == 0)
   
   #check the output if bt is negative
   with pytest.raises(ValueError):
       fn.emittance(wav,-1)
       
   #check if the output length is correct
   assert(len(intensity) == len(wav))
   
# Test for function 'transmittance'

@given(st.floats(0.1,4),st.floats(0.5,1),st.floats(10,15),st.floats(5,9.99),
       st.floats(0,200),st.floats(0,200),st.floats(0.1,1))
@settings(max_examples = 10)

def test_transmittance(start_wav,wav_step,final_wav,wavelength_norm,tau_v,
                       tau_h,ang_mus):
   
   wav, _ = fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm) 
   tau_h = np.random.rand(len(wav))
   tau_v = np.random.rand(len(wav))
   Ts, Th = fn.transmittance(wav,tau_v,tau_h,ang_mus)
   
     
   #check if output is positive
   assert(len(Ts[Ts <=0]) == 0)
   assert(len(Th[Th <=0]) == 0)
   
   #check if length of Th and Ts are equals to length of wav
   assert(len(Ts) == len(wav))
   assert(len(Th) == len(wav))
  
   #check if cosine has value out of the range [-1,1]
   with pytest.raises(ValueError):
        fn.transmittance(wav,tau_v,tau_h,10)
                
# Test for function 'irradiance'
   
@given(st.floats(0.1,4),st.floats(0.5,1),st.floats(10,15),st.floats(5,9.99),
       st.floats(1,8000),st.floats(0,100), st.floats(0.1,1))
@settings(max_examples = 20)

def test_irradiance(start_wav,wav_step,final_wav,wavelength_norm,bt,tau_v,
                    ang_mus):
   
   wav, _ = fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm)
   tau_v = np.random.rand(len(wav))
   E0, E0S, E0N, E0SN = fn.irradiance(wav,wav_step,bt,tau_v,ang_mus)
   
   #check if the E0N and E0SN values are equals or smaller than 1
   assert(len(E0N[E0N <=1]))
   assert(len(E0SN[E0SN <=1]))   
   #check if output is positive
   assert(len(E0[E0 <=0]) == 0)
   assert(len(E0S[E0S <=0]) == 0)
   assert(len(E0N[E0N <=0]) == 0)
   assert(len(E0SN[E0SN <=0]) == 0)
   #check if length of E0,E0S,E0N and E0SN are equals to length of wav
   assert(len(E0) == len(wav))
   assert(len(E0S) == len(wav))
   assert(len(E0N) == len(wav))
   assert(len(E0SN) == len(wav))
   
   #check if cosine has value out of the range [-1,1]
   with pytest.raises(ValueError):
       fn.irradiance(wav,wav_step,bt,tau_v,5)
   #check the output if bt is negative
   with pytest.raises(ValueError):
       fn.irradiance(wav,wav_step,-50,tau_v,ang_mus)
   
   
# Test for function 'ScatteringTot'
   
@given(st.floats(0,100),st.floats(0,100),st.floats(0.1,4),st.floats(0.5,1),
       st.floats(10,15),st.floats(5,9.99))
@settings(max_examples = 20)

def test_ScatteringTot(k_mol_scattering,k_aer_scattering,start_wav,wav_step,
			final_wav,wavelength_norm):
    
    wav, _ = fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm)
    k_mol_scattering = np.random.rand(len(wav))
    k_aer_scattering = np.random.rand(len(wav))
    Stot, ang, Ms, As = fn.ScatteringTot(k_mol_scattering,k_aer_scattering,
                                         wav,wav_step)
    
    #check if output is positive
    assert(len(Stot[Stot <=0]) == 0)
    assert(len(ang[ang < 0]) == 0)
    assert(len(Ms[Ms <=0]) == 0)
    assert(len(As[As <=0]) == 0)

    #check if the output length are correct
    assert(len(Stot) == len(ang))
    assert(len(Ms) == len(ang))
    assert(len(As) == len(ang))  
    
    #check if Stot is composed by zero elements like k_mol/k_aer_scattering 
    assert(len(Stot[Stot ==0]) == len(k_mol_scattering[k_mol_scattering ==0]))
    assert(len(Stot[Stot ==0]) == len(k_aer_scattering[k_aer_scattering ==0]))
    

if __name__ == '__main__':
    pass
   
   
