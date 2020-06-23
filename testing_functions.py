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

@given(st.floats(0,5),st.floats(0,5))
@settings(max_examples = 3)

def test_dir_scatteringTot(dir_mol,dir_aer):
    assert(f.dir_scatteringTot((dir_mol),(dir_aer))>0 and 
           f.dir_scatteringTot((dir_mol),(dir_aer)) <10)
    
    #check that output do not contain negative elements
    assert(len(f.sin.D_tot[0][f.sin.D_tot[0] < 0]) == 0)
    
       
        
# Test for planck function

@given(st.floats(0.3,5),st.floats(1,6000),st.integers(0,100))
@settings(max_examples = 3)

def test_planck(wav,bt,i):
   assert(f.planck(wav,bt,i)>=0)
   with pytest.raises(ValueError):
       f.planck(wav,-1,i)
       
       # check if the output length is correct
       assert(len(f.planck(wav,bt,i)) == i)
       
       #check that output do not contain negative elements
       assert(len(f.sin.intensity[f.sin.intensity < 0]) == 0)

           
# Test for Total scattering function
    
@given(st.floats(0,90),st.floats(0,90))
@settings(max_examples = 3)

def test_ScatteringTot(Mol_scattering,Aer_scattering):
     #check that output do not contain negative elements
    assert(len(f.sin.Stot[f.sin.Stot < 0]) == 0)
    

if __name__ == '__main__':
    pass
    
    
    
    
    