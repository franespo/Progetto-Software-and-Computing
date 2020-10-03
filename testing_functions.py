# Testing functions.py

import numpy as np
import functions as fn
import pytest
from hypothesis import strategies as st
from hypothesis.strategies import tuples
from hypothesis import settings
from hypothesis import given

# Test for wav_matrix
 
@given(st.floats(0,5),st.floats(0,1),st.floats(5,10),st.floats(0,10))
@settings(max_examples = 4)

def test_wav_matrix(start_wav,wav_step,final_wav,wavelength_norm):
    
    # Check if output is positive
    assert(fn.wav_matrix(start_wav,wav_step,final_wav,wavelength_norm)>0)
    
    #check if a ValueError arise if start_wav is smaller than 0
    with pytest.raises(ValueError):
       fn.wav_matrix(-1,0.01,5,2)
        
    #check if ValueError arise if final_wav is smaller than start_wav
    with pytest.raises(ValueError):
       fn.wav_matrix(0.5,0.1,0.4,0.3)
       
       
# Test for radcos

@given(st.floats(0,180),st.floats(0,90))
@settings(max_examples = 4)

def test_radcos(scata,zena):

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
       
@given(st.floats(0,30),st.floats(0,100))
@settings(max_examples = 4)

def test_RayOpticaldepth(Hv,wav):
        
    #check if the output length is the correct
    assert(len(fn.RayOpticaldepth(Hv,wav)) == np.size(wav))
        
    #check if a ValueError arise if Hv is smaller than 0
    with pytest.raises(ValueError):
       fn.RayOpticaldepth(-1,wav)
       
# Test for scatter_coeff
       
@given(st.floats(0,100),st.floats(0,100),st.floats(0,100),
       st.floats(0,100),st.floats(0,50),st.floats(0,20))
@settings(max_examples = 4)

def test_scatter_coeff(tau_vertical,Hv,wav,i_wav,maratio, alpha):
    
    #check if the output length is correct
    assert(len(fn.scatter_coeff(tau_vertical,Hv,wav,i_wav,maratio,
                                alpha) == np.size(wav)))

# Test for Opticaldepth
    
@given(tuples(st.floats(0,100),st.floats(0,100),st.floats(0,100)),
       st.floats(0,100),st.floats(0,100))
@settings(max_examples = 4)

def test_Opticaldepth(k_tot_scattering,Hv,d):
      
    #check if the output length is correct
    assert(len(fn.Opticaldepth(k_tot_scattering,Hv,d) == 
               np.size(k_tot_scattering)))
    
    #check if a ValueError arise if d is 0
    with pytest.raises(ValueError):
        fn.Opticaldepth(k_tot_scattering,Hv,0)    

# Test for dir_scatteringTot

@given(st.floats(0,100),st.floats(0,100),st.floats(0,2*np.pi),st.floats(0.30))
@settings(max_examples = 4)
        
def test_dir_scatteringTot(tau_vertical,Hv,ang_mu,g):
    
    #check if output is positive
    assert(fn.dir_scatteringTot(tau_vertical,Hv,ang_mu,g)>0)
    
    #check if output do not contain negative elements
    assert(len(fn.D_tot[0][fn.D_tot[0] < 0]) == 0)
    
# Test for emittance

@given(st.floats(0,100),st.floats(1,8000))
@settings(max_examples = 4)

def test_emittance(wav,bt):
    
   #check if output is positive
   assert(fn.emittance(wav,bt)>=0)
   
   #testing with a negative value, the value error is detected
   with pytest.raises(ValueError):
       fn.emittance(wav,-1)
       
   #check if the output length is correct
   assert(len(fn.emittance(wav,bt)) == len(wav))
       
   #check if output do not contain negative elements
   assert(len(fn.intensity[fn.intensity < 0]) == 0)
    
# Test for transmittance

@given(st.floats(0,100),st.floats(0,200),st.floats(0,200),st.floats(0,2*np.pi))
@settings(max_examples = 4)

def test_transmittance(wav,tau_v,tau_h,ang_mus):
    
   #check if output is positive
   assert(fn.transmittance(wav,tau_v,tau_h,ang_mus)>=0)
   
   #check if length of th is equal to length of wav
   assert(len(fn.transmittance[fn.th]) == np.size(wav))
   
   #check if length of ts is equal to length of wav
   assert(len(fn.transmittance[fn.ts]) == np.size(wav))
  

# Test for irradiance
   
@given(st.floats(0,100),st.floats(0,10),st.floats(0,8000),st.floats(0,100),
       st.floats(0,2*np.pi))
@settings(max_examples = 4)

def test_irradiance(wav,wav_step,bt,tau_v,ang_mus):
    
   #check if output is positive
   assert(fn.irradiance(wav,wav_step,bt,tau_v,ang_mus)>=0)
    
   #check if length of E0 is equal to length of wav
   assert(len(fn.irradiance[fn.E0]) == np.size(wav))
   
   #check if length of Ts is equal to length of wav
   assert(len(fn.irradiance[fn.Ts]) == np.size(wav))
   
    #check if length of E0S is equal to length of wav
   assert(len(fn.transmittance[fn.E0S]) == np.size(wav))
   

# Test for ScatteringTot
   
@given(st.floats(0,100),st.floats(0,10),st.floats(0,100),st.floats(0,10),
       st.floats(0,20))
@settings(max_examples = 4)

def test_ScatteringTot(k_mol_scattering,k_aer_scattering,wav,wav_step,g): 
    
    #check if output is positive
   assert(fn.ScatteringTot(k_mol_scattering,k_aer_scattering,wav,wav_step,g
                           )>=0)
   
   

   
   
    
    
    
    
    
    
    
    

if __name__ == '__main__':
    pass
    