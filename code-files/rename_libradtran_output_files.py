#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 12:16:24 2024

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

# Directory (in Puhti) containing the .json files (libradtran output files)
json_files_dir = "/fmi/projappl/project_2004400/jamin/data/libradtran_data/output_params/"

# List of all the files in the output_params directory
json_files_list = os.listdir(json_files_dir)

# Loop through all the files
for filename in json_files_list:
    
    # Split the filename (we don't want the file extension to be in the middle of the filename)
    splitted_filename = filename.split(".")
    
    # Create the new filename
    new_filename = f"{splitted_filename[0]}_sza35_phi00_alt0_tau0.15.json"
    
    # Old file path
    old_file_path = os.path.join(json_files_dir, filename)
    
    # New file path
    new_file_path = os.path.join(json_files_dir, new_filename)
    
    # Rename the file
    os.rename(old_file_path, new_file_path)