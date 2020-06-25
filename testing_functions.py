# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 11:13:48 2020

@author: Francesco
"""
import pytest
import functions as f
import hypothesis
from hypothesis import strategies as st
from hypothesis import settings
from hypothesis import given

# Test for total scattering direction function

# Testing on input value of the function with any floats between 0 and 5.
# Three maximum examples
 
@given(st.floats(0,5),st.floats(0,5))
@settings(max_examples = 3)

def test_dir_scatteringTot(dir_mol,dir_aer):
    # Check if output is between 0 and 10
    assert(f.dir_scatteringTot((dir_mol),(dir_aer))>0 and 
           f.dir_scatteringTot((dir_mol),(dir_aer)) <10)
    
    # Check if output do not contain negative elements
    assert(len(f.sin.D_tot[0][f.sin.D_tot[0] < 0]) == 0)
    
       
        
# Test for planck function
    
# Testing on input value of the function with any floats between 0.3 and 5 (for
# wav value), 1 and 6000 (for bt value) and any integers between 0 and 100 (for
# i). Three maximum examples.
 
@given(st.floats(0.3,5),st.floats(1,6000),st.integers(0,100))
@settings(max_examples = 3)

def test_planck(wav,bt,i):
    # Check if output is positive
   assert(f.planck(wav,bt,i)>=0)
   
   # Tesing with a negative value, the value error is detected
   with pytest.raises(ValueError):
       f.planck(wav,-1,i)
       
       # Check if the output length is correct
       assert(len(f.planck(wav,bt,i)) == i)
       
       # Check if output do not contain negative elements
       assert(len(f.sin.intensity[f.sin.intensity < 0]) == 0)

           
# Test for Total scattering function
# Testing on input value of the function with any floats between 0 and 90. 
# Three maximum examples.       
    
@given(st.floats(0,90),st.floats(0,90))
@settings(max_examples = 3)

def test_ScatteringTot(Mol_scattering,Aer_scattering):
    # Check if output is positive
    assert(f.ScatteringTot(Mol_scattering,Aer_scattering)>0)
    
    # Check if output do not contain negative elements
    assert(len(f.sin.Stot[f.sin.Stot < 0]) == 0)
    

if __name__ == '__main__':
    pass
    
    
    
    
    