# #############################################################################
# blaablaa
# ############################################################################# 

# In[]:
# #############################################################################
# Set the working directory
# #############################################################################

import os
os.chdir("/home/jaminkiukkonen/Desktop/work/SummerProject")
os.getcwd()

# In[]:
# #############################################################################
# Packages
# #############################################################################

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import run_libradtran_extended as libis
import configparser
import glob
import json
import pandas as pd

# In[]:
# #############################################################################
# Read in the seasonal temperature profiles
# #############################################################################

def read_netCDF4_data(filename):
    with Dataset(filename) as rootgrp:
        t = rootgrp["data"][:, :]
        return t

# The seasonal temperature profiles
Q1_profiles = read_netCDF4_data("/home/jaminkiukkonen/Desktop/work/SummerProject/data/Q1_profiles.nc")
Q2_profiles = read_netCDF4_data("/home/jaminkiukkonen/Desktop/work/SummerProject/data/Q2_profiles.nc")
Q3_profiles = read_netCDF4_data("/home/jaminkiukkonen/Desktop/work/SummerProject/data/Q3_profiles.nc")
Q4_profiles = read_netCDF4_data("/home/jaminkiukkonen/Desktop/work/SummerProject/data/Q4_profiles.nc")

# print(Q1_profiles.shape)    # (18, 19)
# print(Q2_profiles.shape)    # (18, 19)
# print(Q3_profiles.shape)    # (18, 19)
# print(Q4_profiles.shape)    # (18, 19)

# In[]:
# #############################################################################
# Read in the midlatitude summer file (afglms.dat)
# #############################################################################

# Specify the file path of the file
file_path = "/home/jaminkiukkonen/libRadtran-2.0.5/data/atmmod/afglms.dat"

# Read in the midlatitude summer file as a NumPy array
mls = np.genfromtxt(file_path, dtype="float")

# The 2nd column contains the pressure levels
# print(mls[:, 1].shape)    # (50,)

# The 3rd column contains the temperature profile
# print(mls[:, 2].shape)    # (50,)

# In[]:
    
# print(mls[:, 1] / ((mls[:, 2] * mls[:, 3])))

# In[]:
# #############################################################################
# A function to write over the temperature profile in the midlatitude summer 
# file (afglms.dat). To achieve this, we need to perform interpolation on the 
# seasonal temperature profile so that they match the pressure levels in the
# afglms.dat file. In addition to this, this function recalculates
# the air densities based on the new temperature profile.
# #############################################################################

def modify_mls(T_profile, mls, orig_levels, file_name):
    """
    This function writes over the temperature profile in the afglms.dat file,
    and recalculates the air densities based on the new temperatures.
    
    Parameters:
        T_profile: the new temperature profile that you want to include in the
                   afglms.dat file.
        mls: the afglms.dat file as a NumPy array.
        orig_levels: the original pressure levels (needed for interpolation).
        file_name: the file name that you want to save the modified file as.
        season: the season of the temperature profile (Q1, Q2, Q3 or Q4).
                This is used to name the file so we can trace back to it.
        
    Returns:
        The name of the modified midlatitude summer file.
    """
    new_pressure_levels = mls[:, 1]
    interpolated_temperatures = np.interp(new_pressure_levels,
                                          orig_levels,
                                          T_profile)
    # Change the temperature profile
    mls[:, 2] = interpolated_temperatures
    
    # Recalculate the air densities
    unit_conversion_constant = 10**(-4)
    K = 1.380649 * 10**(-23)     # Boltzmann constant
    p = mls[:, 1]    # Pressure levels (mb)
    T = mls[:, 2]    # Temperatures (K)
    new_ad = (p * unit_conversion_constant) / (K * T)
    # Write over the air density column with the recalculated air densities
    mls[:, 3] = new_ad
    
    # Write the mls array back to a text file
    col_names = " z(km), p(mb), T(K), air(cm-3), o3(cm-3), o2(cm-3), h2o(cm-3), co2(cm-3), no2(cm-3)"
    path = os.getcwd() + "/mls_files/" 
    np.savetxt(path + file_name, 
               mls, 
               header=col_names,
               fmt="%20.6e",
               comments="#")
    return file_name  
    
# In[]:
# #############################################################################
# Generate all the different midlatitude summer files (afglms.dat),
# such that we replace the temperature profile in the afglms.dat file
# with the temperature profile that we created our self.
# We have 18 temperature profiles on each season (quartal),
# which will yield in 4*18=72 different afglms.dat files.
# ############################################################################# 

def create_afglms_files(seasonal_profiles, quartal):
    orig_pressure_levels = np.array([1,5,10,30,50,70,100,125,175,225,300,400,
                                     500,600,700,800,900,950,1000])
    file_names = []
    for k in range(seasonal_profiles.shape[0]):
        T_profile = seasonal_profiles[k, :]
        file_name = "afglms" +  "_" + quartal + "_" + str(k)
        # Call the modify_mls() function to create the file,
        # and save the file name into a list.
        file_names.append(modify_mls(T_profile, 
                                     mls, 
                                     orig_pressure_levels, 
                                     file_name))
    return file_names

# Create the 72 modified afglms (midlatitude summer) files.
# Note that we are saving the file names into a list as well
afglms_file_names = []
afglms_file_names += create_afglms_files(Q1_profiles, "Q1")
afglms_file_names += create_afglms_files(Q2_profiles, "Q2")
afglms_file_names += create_afglms_files(Q3_profiles, "Q3")
afglms_file_names += create_afglms_files(Q4_profiles, "Q4") 