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
import json

# In[]:
# #############################################################################
# Set the working directory
# #############################################################################

os.chdir("/home/jaminkiukkonen/Desktop/work/SummerProject/mls_files/")
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

def change_T_profile(T_profile, mls, orig_levels, file_name):
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
                                     500, 600, 700, 800, 900, 950, 1000])
    file_names = []
    for k in range(seasonal_profiles.shape[0]):
        T_profile = seasonal_profiles[k, :]
        file_name = "afglms" +  "_" + quartal + "_" + str(k)
        # Call the change_T_profile() function to create the file,
        # and while we're at it, save the file name into a list as well.
        file_names.append(change_T_profile(T_profile, 
                                           mls, 
                                           orig_pressure_levels, 
                                           file_name))
    return file_names

# Create the 72 modified afglms (midlatitude summer) files.
# Note that we are saving the file names into a last as well
afglms_file_names = []
afglms_file_names += create_afglms_files(Q1_profiles, "Q1")
afglms_file_names += create_afglms_files(Q2_profiles, "Q2")
afglms_file_names += create_afglms_files(Q3_profiles, "Q3")
afglms_file_names += create_afglms_files(Q4_profiles, "Q4")

# In[]:
# #############################################################################
# Blablabla
# #############################################################################

def run_lib(afglms_file_name: str, deg: float, wlen: float, phi: float, 
            phi0: float, alt: float, sza: float, vza: float, albedo: float, 
            tau: float, zout_type: str, output_quantity_type: str, 
            spectral_reso: str, output_format: int):
    """
    Function info here.
    
    Params:
        afglms_file_name: name of the midlatitude summer file (not full path).
        
        deg: the degree for viewing zenith angle (VZA).
        
        wlen: value for wavelength.
        
        phi: sensor azimuth value.
        
        phi0: azimuth value.
        
        alt: altitude value (distance from the sea level)
        
        sza: solar zenith angle value.
        
        vza: viewing zenith angle.
        
        albedo: the reflectivity of a surface.
        
        tau: the integrated optical thickness can be set to a constant value
        using this parameter.
        
        zout_type: This option is used to specify the output altitudes 
        in km above surface altitude.
        You can also use "toa" for top of atmosphere, 
        "sur" for surface altitude and "cpt" for cold point tropopause.
        
        output_quantity_type: output_quantity in libradtran can be used to
        convert radiances / irradiances to equivalent output quantity,
        where the quantity type can be one of the following:
        brightness, reflectivity, transmittance.
        
        spectral_reso: for both solar and thermal bands, different spectral
        resolutions are available, and can be selected using
        mol_abs_param reptran coarse / medium / fine.
        
        output_format: tells the script the shape of the output string 
        so that the scripts are parsed correctly.
        Shapes are described in libRadtran User's Guide in section 3.1.5.
        
            - output_format = 0 corresponds to string shape were no umu or phi is specified. 
            - output_format = 1 corresponds to string shape were umu is specified but phi is not.
            - output_format = 2 corresponds to string shape were umu and phi are specified.
        
    Returns:
        return info here.
        
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
        
        - aerosol_vulcan and aerosol_haze are set to 1.
        
        - Output altitudes (zout) must be within the range 
        defined in the atmosphere_file.
    """
    # This provides path for the atmosphere_file
    working_directory = os.getcwd()
    path = working_directory + "/" + afglms_file_name
    
    libradtran_data_files_path = "/home/jaminkiukkonen/libRadtran-2.0.5/data"
    libradtran_bin_file_loc = "/home/jaminkiukkonen/libRadtran-2.0.5/bin"
    
    # Convert degrees to radians
    vza = np.deg2rad(deg)
    
    # Input string for the configuration object
    input_str = f"data_files_path {libradtran_data_files_path} \n\
                wavelength {wlen} \n\
                aerosol_default \n\
                atmosphere_file {path} \n\
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
    
    # Save the run configuration into a .ini file
    save_loc_ini = "/home/jaminkiukkonen/Desktop/work/SummerProject/ini_files/"
    with open(save_loc_ini + afglms_file_name + ".ini", "w") as f:
        conf_obj.write(f)
        
    # Save the run config into a JSON file
    def save_config(data, filepath):
        with open(filepath, "w") as json_file:
            json.dump(data, json_file, indent=4)
        
    # Save the parameters of the current run into a dictionary
    param_config = {"atmosphere_file": afglms_file_name,
                    "degree for VZA": deg,
                    "wavelength": wlen,
                    "phi": phi,
                    "phi0": phi0,
                    "altitude": alt,
                    "SZA": sza,
                    "VZA": vza,
                    "albedo": albedo,
                    "tau": tau,
                    "zout": zout_type,
                    "output quantity": output_quantity_type,
                    "spectral resolution": spectral_reso}
    
    # Save the parameter config into a JSON file.
    save_loc_JSON = "/home/jaminkiukkonen/Desktop/work/SummerProject/JSON_files/"
    save_config(param_config, save_loc_JSON + afglms_file_name + ".json")

    
    #results_list = libis.run_conf_files_libradtran(libradtran_bin_file_loc, input_str)
    
    

# In[]:
# #############################################################################
# Test run
# #############################################################################

name = "afglms_Q3_17"

run_lib(name, 45.0, 550.0, 800.0, 180.0, 100.0, 50.0, 50.0, 100.0, 1.0, "toa", "brightness", "coarse", 2)

# In[]:
# #############################################################################
# Now we are ready to run libradtran on the different configs.
# #############################################################################

