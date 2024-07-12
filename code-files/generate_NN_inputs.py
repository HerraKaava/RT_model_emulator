#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 11:21:12 2024

@author: jaminkiukkonen
"""

#%%

import configparser
import re
import pandas as pd
import numpy as np
import glob

#%%

def extract_vals(full_ini_path, key):
    """
    Args:
        full_ini_path -- a full path to the .ini file that you wish to extract values from
        key -- the variable name in the .ini file input string that you wish to extract
    """
    # Initialize the parser
    config = configparser.ConfigParser()
    
    # Read in the .ini file
    config.read(full_ini_path)
    
    # Extract the input string from the main_str section
    input_str = config["main_str"]["input_str"]
    
    if key == "tau":
        tau_match = re.search(r'\baerosol_modify tau set ([\d.]+)', input_str)
        if tau_match:
            return float(tau_match.group(1))
    
    elif key == "atmosphere_file":
        atmosphere_match = re.search(r'\batmosphere_file\s+(\S+)', input_str)
        if atmosphere_match:
            split_path = atmosphere_match.group(1).split("/")
            return split_path[-1]
    
    # Regular expressions to extract the values of interest
    match = re.search(rf'\b{key}\b\s+([\d.]+)', input_str)
    
    if match:
        return float(match.group(1))
    else:
        raise ValueError(f"Invalid key: {key} not found")

#%%

def construct_matrix(full_ini_path):
    """
    Constructs an input matrix for the NN by extracting the following values
    from the .ini files containing the run configs: phi, phi0, sza, altitude, tau, umu.
    The seasonal atmospheric profile is used as an identifier,
    and its values with their corresponding pressure levels are concatenated to the input matrix as well.
    
    Args:
        full_ini_path -- a full path to the .ini file that you wish to extract values from
        
    Returns:
        input_matrix -- an input matrix for the neural network.
    """
    mls_folder_Path = "/home/jaminkiukkonen/Desktop/tensorflow/mls_files/"
    
    phi = extract_vals(full_ini_path, key="phi")
    phi0 = extract_vals(full_ini_path, key="phi0")
    sza = extract_vals(full_ini_path, key="sza")
    altitude = extract_vals(full_ini_path, key="altitude")
    tau = extract_vals(full_ini_path, key="tau")
    umu = extract_vals(full_ini_path, key="umu")
    atmosphere_file = extract_vals(full_ini_path, key="atmosphere_file")
    
    # Read in the midlatitude summer file as a NumPy array
    mls_arr = np.genfromtxt(mls_folder_Path + atmosphere_file, dtype="float")
    
    # The third column contains the temperature profile
    T_profile = mls_arr[:, 2]
    
    # The second column contains the pressure levels (these are used as column names)
    pressure_levels = [int(lvl) for lvl in mls_arr[:, 1]]
    
    # Create a dataframe (row vector) of temperatures where the pressure levels are the column names
    df_T = pd.DataFrame([T_profile], columns=[str(level) for level in pressure_levels])
    
    # Create a dataframe (row vector) of the variables of interest
    df_vars = pd.DataFrame({"atmosphere_file": atmosphere_file,
                            "phi": [phi],
                            "phi0": [phi0],
                            "sza": [sza],
                            "altitude": [altitude],
                            "tau": [tau],
                            "umu": [umu]})
    
    # Concatenate the two dataframes to form a one bigger dataframe (row vector)
    input_matrix = pd.concat([df_vars, df_T], axis=1)
    
    return input_matrix

#%%

# Path to the folder containing the .ini files
ini_folder_path = "/home/jaminkiukkonen/Desktop/tensorflow/ini_files/"

# These files contain the unique run configurations that we need
filtered_ini_paths = sorted(glob.glob(ini_folder_path + "*conf_obj1*"))

#%%

# Initialize an empty DataFrame to store the concatenated results
X = pd.DataFrame()

# Loop over each file path, call construct_matrix, and concatenate the results
for ini_path in filtered_ini_paths:
    M = construct_matrix(ini_path)
    X = pd.concat([X, M], axis=0, ignore_index=True)
    
#%%

    
    

    
    
    
    
    





