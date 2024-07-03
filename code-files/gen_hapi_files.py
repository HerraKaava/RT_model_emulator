#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 09:02:26 2024

@author: jaminkiukkonen
"""

# In[]:
# #############################################################################
# Packages
# #############################################################################

import numpy as np
import shutil
import os

# In[]:
    
os.getcwd()

# In[]:
# #############################################################################
# Define paths
# #############################################################################
    
# Path to original midlatitude summer file
path_original_mls = "/home/jaminkiukkonen/libRadtran-2.0.5/data/atmmod/afglms.dat"

# Directory containing the modified midlatitude summer files
dir_modified_mls = "/home/jaminkiukkonen/Desktop/work/SummerProject/mls_files/"

# hapi_tests_dev.py script path
hapi_dev_path = "/home/jaminkiukkonen/Desktop/run-libradtran/hapi_tests_dev.py"

# In[]:   

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

# In[]:
# #############################################################################
# Generate files
# #############################################################################
    
for i in range(len(afglms_file_names)):
    
    # Get the file path for the ith modified mls file
    source_file = dir_modified_mls + afglms_file_names[i]
    
    # Overwrite the contents of the original mls with the modified mls
    shutil.copyfile(source_file, path_original_mls)
    
    # Output file name provided by hapi_tests_dev.py
    output_name = "uvspec_afglms_libis_test_500-800nm.nc"
    
    # Define the new output file name with an identifier
    new_output_name = f"uvspec_{afglms_file_names[i]}.nc"
    
    # Directory to save the output files
    save_dir = "/home/jaminkiukkonen/Desktop/work/SummerProject/hapi_output_files/"
    
    # Run the script and move the output file to desired location.
    # Using os.system() allows you to execute scripts 
    # as if you were running it from the terminal.
    os.system(f"python3 {hapi_dev_path} && mv {output_name} {save_dir}/{new_output_name}")







