#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 11:37:22 2024

@author: jaminkiukkonen
"""

#%%

import os

#%%
# #############################################################################
# Renaming the files that got generated in Puhti with the first run.
# The files will be modified such that identifiers for the files will be
# sza=35, phi0=0, altitude=0, and tau=0.15. In addition to this, the files
# already contain the seasonal atmospheric file name identifier.
# #############################################################################

# Directory containing the .ini files
ini_file_dir = "/home/jaminkiukkonen/Desktop/puhti_data/ini_files/"

# List of all the files in the ini_files directory
ini_files_list = os.listdir(ini_file_dir)

# Loop through all the files
for filename in ini_files_list:
    
    # Split the filename (we don't want the file extension to be in the middle of the filename)
    splitted_filename = filename.split(".")
    
    # Create the new filename
    new_filename = f"{splitted_filename[0]}_sza35_vza0_phi0_phi00_alt0_tau0.15.ini"
    
    # Old file path
    old_file_path = os.path.join(ini_file_dir, filename)
    
    # New file path
    new_file_path = os.path.join(ini_file_dir, new_filename)
    
    # Rename the file
    os.rename(old_file_path, new_file_path)
    
    
    
    
    
    
    
    
    
    