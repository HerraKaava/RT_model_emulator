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
# afglms.dat file.
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
        The file path of the modified midlatitude summer file.
    """
    new_pressure_levels = mls[:, 1]
    interpolated_temperatures = np.interp(new_pressure_levels,
                                          orig_levels,
                                          T_profile)
    # Change the temperature profile in mls
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
                                     500,600,700,800,900,950,1000])
    file_names = []
    for k in range(seasonal_profiles.shape[0]):
        T_profile = seasonal_profiles[k, :]
        file_name = "afglms" +  "_" + quartal + "_" + str(k)
        # Call the change_T_profile() function to create the file,
        # and while we're at it, save the file name into a list as well.
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

# In[]:
# #############################################################################
# Blablabla
# #############################################################################

def run_lib(afglms_file_name: str, vza_deg: float, wlen_min: float, wlen_max: float, 
            phi: float, phi0: float, alt: float, sza: float, vza: float,
            albedo1: float, albedo2: float,tau: float, output_format: int, spectral_reso: str):
    """
    Function info here.
    
    Params:
        afglms_file_name: name of the midlatitude summer file (not full path).
        
        vza_deg: the degree for viewing zenith angle (VZA).
        
        wlen: value for wavelength.
        
        phi: sensor azimuth angle
        
        phi0: solar azimuth angle
        
        alt: altitude value (distance from the sea level)
        
        sza: solar zenith angle value.
        
        vza: viewing zenith angle.
        
        albedo1: the reflectivity of a surface.
        
        albedo2: the reflectivity of a surface.
        
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
    albedo0 = 0.0
    
    # This provides path for the atmosphere_file
    path_afglms = os.getcwd() + "/mls_files/" + afglms_file_name
    
    libradtran_data_files_path = "/home/jaminkiukkonen/libRadtran-2.0.5/data"
    libradtran_bin_file_loc = "/home/jaminkiukkonen/libRadtran-2.0.5/bin"
    
    # Convert degrees to radians
    vza = np.deg2rad(vza_deg)
    
    def create_input_str(albedo, output_quantity_type, zout_type):
        # Input string for the configuration object
        high_reso_path = "/home/jaminkiukkonen/Desktop/work/SummerProject/uvspec_afglms_libis_flex_500-800nm.nc"
        input_str = f"data_files_path {libradtran_data_files_path} \n\
                    wavelength {wlen_min} {wlen_max} \n\
                    aerosol_default \n\
                    atmosphere_file {path_afglms} \n\
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
                    #mol_tau_file abs {high_reso_path} \n\
                    mol_abs_param reptran {spectral_reso} \n\
                    rte_solver disort \n\
                    quiet"                
        if output_quantity_type != "radiance":
            input_str += f"\noutput_quantity {output_quantity_type}" 
        return input_str
    
    input_str1 = create_input_str(albedo1, "radiance", "toa")
    input_str2 = create_input_str(albedo1, "radiance", "sur")
    input_str3 = create_input_str(albedo1, "reflectivity", "sur")
    input_str4 = create_input_str(albedo0, "reflectivity", "toa")
    input_str5 = create_input_str(albedo1, "reflectivity", "toa")
    input_str6 = create_input_str(albedo2, "reflectivity", "toa")
    input_str7 = create_input_str(albedo0, "radiance", "toa")
    
    def create_conf_obj(input_params):
        conf_obj = configparser.ConfigParser()
        conf_obj["parser_options"] = {"output_format": output_format,
                                      "umu_length": vza.size}
        conf_obj["main_str"] = {"input_str": input_params}
        return conf_obj
        
    conf_obj1 = create_conf_obj(input_str1)
    conf_obj2 = create_conf_obj(input_str2)
    conf_obj3 = create_conf_obj(input_str3)
    conf_obj4 = create_conf_obj(input_str4)
    conf_obj5 = create_conf_obj(input_str5)
    conf_obj6 = create_conf_obj(input_str6)
    conf_obj7 = create_conf_obj(input_str7)
    
    conf_obj_list = [conf_obj1, conf_obj2, conf_obj3, conf_obj4, conf_obj5, conf_obj6, conf_obj7]
    
    for idx, conf_obj in enumerate(conf_obj_list):
        # Save the run configuration into a .ini file
        path_ini = os.getcwd() + "/ini_files/"
        with open(path_ini + afglms_file_name + "_conf_obj" + str(idx + 1) + ".ini", "w") as f:
            conf_obj.write(f)
            
    # Selecting the .ini files
    ini_files_path = "/home/jaminkiukkonen/Desktop/work/SummerProject/ini_files/"
    pattern = "*.ini"
    input_files_list = sorted(glob.glob(ini_files_path + pattern))
    
    # Libradtran outputs
    results_list = libis.run_conf_files_libradtran(libradtran_bin_file_loc, input_files_list)

    def save_output_params(results):
        # Dictionary to store the output params
        data = {
            "Wavelength": [],
            "TOA_RAD": [],
            "Edif0": [],
            "Edir0": [],
            "Tdir": [],
            "Tdif": [],
            "Rho0": [],
            "Rho1": [],
            "Rho2": [],
            "Spherical_albedo": [],
            "albedo1": np.repeat(albedo1, repeats=len(results[0])),
            "albedo2": np.repeat(albedo2, repeats=len(results[0])),
            "Path_rad": []
        }
        ########## uu (TOA_RAD) ##########
        # Extract lambda (wavelength) and uu (TOA_RAD) from the 1st list
        for d in results[0]:
            lam = d["lambda"]
            toa_rad = d["uu(umu,phi)"][0][0]
            data["Wavelength"].append(lam)
            data["TOA_RAD"].append(toa_rad)

        ########## edn (edif0), edir (edir0) ##########
        # Extract edn and edir from the 2nd list
        for d in results[1]:
            edif0 = d["edn"]
            edir0 = d["edir"]
            data["Edif0"].append(edif0)
            data["Edir0"].append(edir0)
            
        ########## edir (tdir), edn (tdif) ##########
        # Extract edir and edn from the 3rd list
        for d in results[2]:
            tdir = d["edir"]
            tdif = d["edn"]
            data["Tdir"].append(tdir)
            data["Tdif"].append(tdif)
        
        ########## Rho0, Rho1, Rho2 ##########
        # Extract rho0 from the 4th list
        for d in results[3]:
            rho0 = d["uu(umu,phi)"][0][0]
            data["Rho0"].append(rho0)
            
        # Extract rho1 from the 5th list
        for d in results[4]:
            rho1 = d["uu(umu,phi)"][0][0]
            data["Rho1"].append(rho1)
            
        # Extract rho2 from the 6th list
        for d in results[5]:
            rho2 = d["uu(umu,phi)"][0][0]
            data["Rho2"].append(rho2)
            
        ########## Spherical albedo ##########
        numerator1 = albedo1 * (np.array(data["Rho2"]) - np.array(data["Rho0"]))
        numerator2 = albedo2 * (np.array(data["Rho1"]) - np.array(data["Rho0"]))
        denominator = albedo2 * albedo1 * (np.array(data["Rho2"]) - np.array(data["Rho1"]))
        spherical_albedo = (numerator1 - numerator2) / denominator
        data["Spherical_albedo"] = spherical_albedo
            
        ########## Path radiance ##########
        # Extract path radiance from the 7th (last) list
        for d in results[6]:
            path_rad = d["uu(umu,phi)"][0][0]
            data["Path_rad"].append(path_rad)
        
        # Set display precision for floating-point numbers
        pd.set_option('display.precision', 10)
        
        # Create a dataframe
        df = pd.DataFrame(data)
        
        # Save df DataFrame into JSON file
        output_params_path = "/home/jaminkiukkonen/Desktop/work/SummerProject/output_params/"
        df.to_json(output_params_path + afglms_file_name + ".json", orient="columns", indent=4)
    
    save_output_params(results_list)

# In[]:
# #############################################################################
# Test run
# #############################################################################

name = "afglms_Q3_17"

y = run_lib(afglms_file_name=name, 
            vza_deg=45.0, 
            wlen_min=500.0, 
            wlen_max=800.0, 
            phi=180.0, 
            phi0=180.0, 
            alt=1.0, 
            sza=40.0, 
            vza=15.0, 
            albedo1=0.15,
            albedo2=0.25, 
            tau=1.0,
            output_format=2, 
            spectral_reso="coarse")
