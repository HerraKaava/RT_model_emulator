#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 11:32:59 2024

@author: jaminkiukkonen
"""

from fixed_lib_run import run_lib
import sys

#%%

def run_RT(start_idx):
    """
    This function runs the radiative transfer (RT) model.
    
    Params:
        start_idx: starting index of the loop.
        end_idx: ending index of the loop.
        
    Notes:
        - The starting and ending indices are needed for the loop because
        one run of the RT model takes several hours, so the runs are 
        parallelized in Puhti such that we launch several runs simultaneously
        from terminal. In other words, start_idx and end_idx controls which
        atmospheric files and hapi files are fed into libradtran as inputs.
    """
    afglms_file_names = [
         'afglms_Q1_0',
         'afglms_Q1_1',
         'afglms_Q1_2',
         'afglms_Q1_3',
         'afglms_Q1_4',
         'afglms_Q1_5',
         'afglms_Q1_6',
         'afglms_Q1_7',
         'afglms_Q1_8',
         'afglms_Q1_9',
         'afglms_Q1_10',
         'afglms_Q1_11',
         'afglms_Q1_12',
         'afglms_Q1_13',
         'afglms_Q1_14',
         'afglms_Q1_15',
         'afglms_Q1_16',
         'afglms_Q1_17',
         'afglms_Q2_0',
         'afglms_Q2_1',
         'afglms_Q2_2',
         'afglms_Q2_3',
         'afglms_Q2_4',
         'afglms_Q2_5',
         'afglms_Q2_6',
         'afglms_Q2_7',
         'afglms_Q2_8',
         'afglms_Q2_9',
         'afglms_Q2_10',
         'afglms_Q2_11',
         'afglms_Q2_12',
         'afglms_Q2_13',
         'afglms_Q2_14',
         'afglms_Q2_15',
         'afglms_Q2_16',
         'afglms_Q2_17',
         'afglms_Q3_0',
         'afglms_Q3_1',
         'afglms_Q3_2',
         'afglms_Q3_3',
         'afglms_Q3_4',
         'afglms_Q3_5',
         'afglms_Q3_6',
         'afglms_Q3_7',
         'afglms_Q3_8',
         'afglms_Q3_9',
         'afglms_Q3_10',
         'afglms_Q3_11',
         'afglms_Q3_12',
         'afglms_Q3_13',
         'afglms_Q3_14',
         'afglms_Q3_15',
         'afglms_Q3_16',
         'afglms_Q3_17',
         'afglms_Q4_0',
         'afglms_Q4_1',
         'afglms_Q4_2',
         'afglms_Q4_3',
         'afglms_Q4_4',
         'afglms_Q4_5',
         'afglms_Q4_6',
         'afglms_Q4_7',
         'afglms_Q4_8',
         'afglms_Q4_9',
         'afglms_Q4_10',
         'afglms_Q4_11',
         'afglms_Q4_12',
         'afglms_Q4_13',
         'afglms_Q4_14',
         'afglms_Q4_15',
         'afglms_Q4_16',
         'afglms_Q4_17'
     ]
    
    wlen_min, wlen_max = 500.0008435, 799.9996703
    sza = 35
    vza = 0
    phi = 0
    phi0 = 0
    alt = 0
    albedo1 = 0.15
    albedo2 = 0.25
    tau = 0.15
                
    # For Puhti
    start_idx = int(start_idx) - 1
    end_idx = start_idx + 1
    
    current_slice = afglms_file_names[start_idx:end_idx]
    for file in current_slice:
        run_lib(file, 
                wlen_min, 
                wlen_max, 
                phi, 
                phi0, 
                alt, 
                sza, 
                vza, 
                albedo1, 
                albedo2, 
                tau)


#%%

if __name__ == "__main__":
    
    run_RT(sys.argv[1])
    
    
    
    
