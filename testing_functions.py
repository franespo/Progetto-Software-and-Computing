# Testing functions.py

import numpy as np
import functions as fn
import pytest
from hypothesis import strategies as st
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
       
       
@given(st.floats(0,180),st.floats(0,90))
@settings(max_examples = 4)

# Test for radcos

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
def test_RayOpticaldepth(Hv,wav):
    
    #check if the output length is the correct
    assert(len(fn.RayOpticaldepth(Hv,wav)) == np.size(wav))
        
    #check if a ValueError arise if Hv is smaller than 0
    with pytest.raises(ValueError):
       fn.RayOpticaldepth(-1,wav)
       









if __name__ == '__main__':
    pass   