# #############################################################################
# blaablaa
# ############################################################################# 

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
import os

# In[]:
# #############################################################################
# Set the working directory
# #############################################################################

os.chdir("/home/jaminkiukkonen/Desktop/run-libradtran/mls_files/")
os.getcwd()

# In[]:
# #############################################################################
# Read in the seasonal temperature profiles
# #############################################################################

def read_netCDF4_data(filename):
    with Dataset(filename) as rootgrp:
        t = rootgrp["data"][:, :]
        return t

# A seasonal temperature profile
Q1_profiles = read_netCDF4_data("/home/jaminkiukkonen/Desktop/work/SummerProject/data/Q1_profiles.nc")
Q2_profiles = read_netCDF4_data("/home/jaminkiukkonen/Desktop/work/SummerProject/data/Q2_profiles.nc")
Q3_profiles = read_netCDF4_data("/home/jaminkiukkonen/Desktop/work/SummerProject/data/Q3_profiles.nc")
Q4_profiles = read_netCDF4_data("/home/jaminkiukkonen/Desktop/work/SummerProject/data/Q4_profiles.nc")

print(Q1_profiles.shape)    # (18, 19)
print(Q2_profiles.shape)    # (18, 19)
print(Q3_profiles.shape)    # (18, 19)
print(Q4_profiles.shape)    # (18, 19)

# In[]:
# #############################################################################
# Read in the midlatitude summer file (afglms.dat)
# #############################################################################

# Specify the file path of the file
file_path = "/home/jaminkiukkonen/libRadtran-2.0.5/data/atmmod/afglms.dat"

# Read in the midlatitude summer file as a NumPy array
mls = np.genfromtxt(file_path, dtype="float")

# The 2nd column contains the pressure levels
print(mls[:, 1].shape)    # (50,)

# The 3rd column contains the temperature profile
print(mls[:, 2].shape)    # (50,)

# In[]:
# #############################################################################
# A function to write over the temperature profile in the midlatitude summer 
# file (afglms.dat). To achieve this, we need to perform interpolation on the 
# seasonal temperature profile so that they match the pressure levels in the
# afglms.dat file.
# #############################################################################

def change_T_profile(T_profile, mls, orig_levels, file_name, season):
    """
    This function writes over the temperature profile in the afglms.dat file.
    
    Parameters:
        T_profile: the new temperature profile that you want to include in the
                   afglms.dat file.
        mls: the afglms.dat file as a NumPy array.
        orig_levels: the original pressure levels (needed for interpolation).
        file_name: the file name that you want to save the modified file as.
        season: the season of the temperature profile (Q1, Q2, Q3 or Q4).
                This is used to name the file so we can trace back to it.
        
    Returns:
        The file path of the modified midlatitude summer file.
    """
    new_pressure_levels = mls[:, 1]
    interpolated_temperatures = np.interp(new_pressure_levels,
                                          orig_levels,
                                          T_profile)
    # Change the temperature profile in mls
    mls[:, 2] = interpolated_temperatures
    
    # Write the mls array back to a text file
    col_names = " z(km), p(mb), T(K), air(cm-3), o3(cm-3), o2(cm-3), h2o(cm-3), co2(cm-3), no2(cm-3)"
    np.savetxt(file_name, 
               mls, 
               header=col_names, 
               delimiter=",",
               fmt="%20.6e",
               comments="#")
    modified_afglms = "/home/jaminkiukkonen/libRadtran-2.0.5/data/atmmod/modified_afglms.dat"
    return modified_afglms

# In[]:
# #############################################################################
# Demo on how to generate the files on a loop
# ############################################################################# 
    
orig_levels = np.array([1, 5, 10, 30, 50, 70, 100, 125, 175, 225, 300, 
                        400, 500, 600, 700, 800, 900, 950, 1000])

for k in range(Q1_profiles.shape[0]):
    T_profile = Q1_profiles[k, :]
    season = "Q1"
    file_name = "afglms" +  "_" + season + "_" + str(k)
    change_T_profile(T_profile, mls, orig_levels, file_name, season)
    
# In[]:
# #############################################################################
# Blablabla
# #############################################################################

def run_lib(afglms_file, deg, wlen, phi, phi0, alt, sza, vza, albedo, 
            zout_type, tau, output_quantity_type, spectral_reso, 
            output_format, ini_file_name, netcdf4_file_name):
    """
    Function info here.
    
    Params:
        deg: the degree for viewing zenith angle (VZA).
        wlen: value for wavelength.
        phi: sensor azimuth value.
        phi0: azimuth value.
        alt: altitude value (distance from the sea level)
        sza: solar zenith angle value.
        vza: viewing zenith angle.
        albedo: the reflectivity of a surface.
        zout_type: zout gives the sensor altitude altitude above the ground.
                   possible options for zout are toa (top of atmosphere),
                   and ... ??
        tau: the integrated optical thickness can be set to a constant value
        using this parameter.
        output_quantity_type: output_quantity in libradtran can be used to
        convert radiances / irradiances to equivalent output quantity,
        where the quantity type can be one of the following:
        brightness, reflectivity, transmittance.
        spectral_reso: for both solar and thermal bands, different spectral
        resolutions are available, and can be selected using
        mol_abs_param reptran coarse / medium / fine.
        This parameter selects one of those (coarse, medium, fine).
        output_format: needs to be one of the following: 0, 1, 2.
        ini_file_name: name for the .ini file for saving the run config.
        netcdf4_file_name: name for the netCDF4 file for saving the run config.
        
    Returns:
        
    Notes:
        - VZA is converted into radians inside the function.
        - phi=phi0 indicates that the sensor looks into the direction of the sun.
        - phi-phi0=180(deg) means that the sun is in the back of the sensor.
        - Altitude specifies the height at which the atmospheric parameters,
        such as pressure and temperature, are defined for the radiative
        transfer calculations.
        - The azimuth (phi0) is only required for radiance calculations.
        - Viewing direction (umu, phi) and sun position (sza, phi0)
        are defined with respect to the sensor position specified by zout.
        - output_format = 0 corresponds to string shape were no umu or phi is specified. 
        - output_format = 1 corresponds to string shape were umu is specified but phi is not.
        - output_format = 2 corresponds to string shape were umu and phi are specified.
        - output_format tells the script the shape of the output string so that 
        scripts are parsed correctly. 
        - Shapes are described in libRadtran User's Guide in section 3.1.5.
        - aerosol_vulcan and aerosol_haze are set to 1.
        
    """
    libradtran_data_files_path = "/home/jaminkiukkonen/libRadtran-2.0.5/data"
    libradtran_bin_file_loc = "/home/jaminkiukkonen/libRadtran-2.0.5/bin"
    
    # Convert degrees to radians
    vza = np.deg2rad(deg)
    
    # Input string for the configuration object
    input_str = f"data_files_path {libradtran_data_files_path} \n\
                wavelength {wlen} \n\
                aerosol_default \n\
                atmosphere_file {afglms_file} \n\
                aerosol_vulcan 1 \n\
                aerosol_haze 1 \n\
                phi {phi} \n\
                phi0 {phi0} \n\
                altitude {alt} \n\
                sza {sza} \n\
                umu {vza} \n\
                albedo {albedo} \n\
                zout {zout_type} \n\
                aerosol_modify tau set {tau} \n\
                output_quantity {output_quantity_type} \n\
                mol_abs_param reptran {spectral_reso} \n\
                rte_solver disort \n\
                verbose"
    
    conf_obj = configparser.ConfigParser()
    conf_obj["parser_options"] = {"output_format": output_format,
                                "umu_length": vza.size}
    conf_obj["main_str"] = {"input_str": input_str}
    
    # Saving the run configuration into .ini format
    with open(ini_file_name, "W") as f:
        conf_obj.write(f)
        
    # Saving the run config into netCDF4 file
    def save_config(data, filename):
        with Dataset(f"{filename}.nc", "w", format="NETCDF4") as nc:
            # Description
            nc.description = "The values of the data file correspond to these variables in this order:\
            afglms_file, deg, wlen, phi, phi0, alt, sza, vza, albedo,\
            zout_type, tau, output_quantity_type, spectral_reso" 
            # Dimensions
            nc.createDimension("columns", len(data))
            # Variables
            nc_vars = nc.createVariable("params", "f8", ("columns",))
            # Assign the values
            nc_vars[:] = data
            
    results_list = libis.run_conf_files_libradtran(libradtran_bin_file_loc, input_str)
    
    # WORKING PROGESS, SCRIPT NOT READY! BACKING UP IN GITHUB.
    

# In[]:
# #############################################################################
# Test run
# #############################################################################

