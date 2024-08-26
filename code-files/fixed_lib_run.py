#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 13:14:52 2024

@author: jaminkiukkonen
"""

#%%

import numpy as np
import run_libradtran_extended as libis
import configparser
import pandas as pd

#%%

def run_lib(afglms_file_name: str, wlen_min: float, wlen_max: float, 
            phi: float, phi0: float, alt: float, sza: float, vza: float,
            albedo1: float, albedo2: float, tau: float):
    """
    Function info here.
    
    Params:
        afglms_file_name: name of the midlatitude summer file (not full path).
        
        vza_deg: the degree for viewing zenith angle (VZA).
        
        wlen_min: minimum value for wavelength.
        
        wlen_max: maximum value for wavelength.
        
        phi: sensor azimuth angle.
        
        phi0: solar azimuth angle.
        
        alt: altitude value (distance from the sea level).
        
        sza: solar zenith angle value.
        
        vza: viewing zenith angle.
        
        albedo1: the reflectivity of a surface.
        
        albedo2: the reflectivity of a surface.
        
        tau: the integrated optical thickness can be set to a constant value
        using this parameter.
        
    Returns:
        None.
        
    Notes:
        - phi=phi0 indicates that the sensor looks into the direction of the sun.
        
        - phi-phi0=180(deg) means that the sun is in the back of the sensor.
        
        - Altitude specifies the height at which the atmospheric parameters,
        such as pressure and temperature, are defined for the radiative
        transfer calculations.
        
        - The azimuth (phi0) is only required for radiance calculations.
        
        - Viewing direction (umu, phi) and sun position (sza, phi0)
        are defined with respect to the sensor position specified by zout.
        
        - Output altitudes (zout) must be within the range 
        defined in the atmosphere_file.
        
        - output_quantity_type param in libradtran can be used to
        convert radiances / irradiances to equivalent output quantity,
        where the quantity type can be one of the following:
        brightness, reflectivity, transmittance.
        
        - zout_type param is used to specify the output altitudes 
        in km above surface altitude.
        You can also use "toa" for top of atmosphere, 
        "sur" for surface altitude and "cpt" for cold point tropopause.
    """
    albedo0 = 0
    
    # This provides path for the atmosphere_file
    path_afglms = "/fmi/projappl/project_2001985/jamin/data/libradtran_data/mls_files/" + afglms_file_name
    
    libradtran_general_path = "/fmi/projappl/project_2001985/APPS/libRadtran/libRadtran-2.0.4"
    libradtran_data_files_path = f"{libradtran_general_path}/data"
    libradtran_bin_file_loc = f"{libradtran_general_path}/bin"
    libradtran_kurudz_full_path = f"{libradtran_general_path}/data/solar_flux/kurudz_full.dat"
    
    # High resolution files path
    hapi_file_path = f"/fmi/projappl/project_2001985/jamin/data/libradtran_data/hapi_output_files/uvspec_{afglms_file_name}.nc"
    
    # ini file save location
    ini_files_save_loc = "/fmi/projappl/project_2001985/jamin/data/libradtran_data/ini_files/"
    
    # Saving location of libradtran output
    output_file_loc_JSON = "/fmi/projappl/project_2001985/jamin/data/libradtran_data/output_params/"
    
    # Convert degrees to radians
    vza_deg = vza
    vza = np.cos(np.deg2rad(vza))
    sza_radian = np.deg2rad(sza)
    
    def create_input_str(albedo, output_quantity_type, zout_type):
        # Input string for the configuration object
        input_str = f"data_files_path {libradtran_data_files_path} \n\
                    wavelength {wlen_min} {wlen_max} \n\
                    aerosol_default \n\
                    source solar {libradtran_kurudz_full_path} \n\
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
                    mol_tau_file abs {hapi_file_path} \n\
                    rte_solver disort \n\
                    quiet"                
        if output_quantity_type != "radiance":
            input_str += f"\noutput_quantity {output_quantity_type}" 
        return input_str
    
    def up_input_str(albedo, output_quantity_type, zout_type):
        # Input string for the configuration object
        input_str = f"data_files_path {libradtran_data_files_path} \n\
                    wavelength {wlen_min} {wlen_max} \n\
                    aerosol_default \n\
                    source solar {libradtran_kurudz_full_path} \n\
                    atmosphere_file {path_afglms} \n\
                    aerosol_vulcan 1 \n\
                    aerosol_haze 1 \n\
                    phi0 {phi0} \n\
                    altitude {alt} \n\
                    sza {vza_deg} \n\
                    albedo {albedo} \n\
                    zout {zout_type} \n\
                    aerosol_modify tau set {tau} \n\
                    mol_tau_file abs {hapi_file_path} \n\
                    rte_solver disort \n\
                    reverse_atmosphere \n\
                    quiet"
        if output_quantity_type != "radiance":
            input_str += f"\noutput_quantity {output_quantity_type}"
        return input_str

    def down_input_str(albedo, output_quantity_type, zout_type):
        # Input string for the configuration object
        input_str = f"data_files_path {libradtran_data_files_path} \n\
                    wavelength {wlen_min} {wlen_max} \n\
                    aerosol_default \n\
                    source solar {libradtran_kurudz_full_path} \n\
                    atmosphere_file {path_afglms} \n\
                    aerosol_vulcan 1 \n\
                    aerosol_haze 1 \n\
                    phi0 {phi0} \n\
                    altitude {alt} \n\
                    sza {sza} \n\
                    albedo {albedo} \n\
                    zout {zout_type} \n\
                    aerosol_modify tau set {tau} \n\
                    mol_tau_file abs {hapi_file_path} \n\
                    rte_solver disort \n\
                    quiet"
        if output_quantity_type != "radiance":
            input_str += f"\noutput_quantity {output_quantity_type}"
        return input_str
    
    interrogation_conf_0 = create_input_str(albedo0, "reflectivity", "toa")
    interrogation_conf_1 = create_input_str(albedo1, "reflectivity", "toa")
    interrogation_conf_2 = create_input_str(albedo2, "reflectivity", "toa")
    down_conf = down_input_str(albedo0, "reflectivity", "sur")
    up_conf = up_input_str(albedo0, "reflectivity", "sur")
      
    def create_conf_obj(input_params, output_format, vza_len):
        conf_obj = configparser.ConfigParser()
        conf_obj["parser_options"] = {"output_format": output_format,
                                      "umu_length": vza_len}
        conf_obj["main_str"] = {"input_str": input_params}
        return conf_obj
        
    conf_obj1 = create_conf_obj(interrogation_conf_0, 2, 1)
    conf_obj2 = create_conf_obj(interrogation_conf_1, 2, 1)
    conf_obj3 = create_conf_obj(interrogation_conf_2, 2, 1)
    conf_obj4 = create_conf_obj(down_conf, 0, 1)
    conf_obj5 = create_conf_obj(up_conf, 0, 1)
    
    conf_obj_list = [conf_obj1, conf_obj2, conf_obj3, conf_obj4, conf_obj5]
    
    for idx, conf_obj in enumerate(conf_obj_list):
        ini_file_name = f"{afglms_file_name}_conf_obj{str(idx+1)}_sza{sza}_vza{vza}_phi{phi}_phi0{phi0}_alt{alt}_tau{tau}.ini"
        # Save the run configuration into a .ini file
        with open(ini_files_save_loc + ini_file_name, "w") as f:
            conf_obj.write(f)
    
    # Libradtran outputs
    results_list = libis.run_conf_files_libradtran(libradtran_bin_file_loc, conf_obj_list)

    def save_output_params(results):
        # Dictionary to store the output params
        data = {
            "wavelength": [],
            "rho0": [],
            "rho1": [],
            "rho2": [],
            "tdir_down": [],
            "tdif_down": [],
            "tdir_up": [],
            "tdif_up": [],
            "spherical_albedo": [],
            "edir": [],
            "edif": [],
            "path_rad": [],
            "albedo1": np.repeat(albedo1, repeats=len(results[0])),
            "albedo2": np.repeat(albedo2, repeats=len(results[0]))
        }
        ########## lambda (wavelength), uu (Rho0) ##########
        for d in results[0]:
            lam = d["lambda"]
            Rho0 = d["uu(umu,phi)"][0][0]
            data["wavelength"].append(lam)
            data["rho0"].append(Rho0)

        ########## uu (Rho1) ##########
        for d in results[1]:
            Rho1 = d["uu(umu,phi)"][0][0]
            data["rho1"].append(Rho1)

        ########## uu(Rho2) ##########
        for d in results[2]:
            Rho2 = d["uu(umu,phi)"][0][0]
            data["rho2"].append(Rho2)

        ########## down ##########
        for d in results[3]:
            tdir_down = d["edir"]
            tdif_down = d["edn"]
            data["tdir_down"].append(tdir_down)
            data["tdif_down"].append(tdif_down)

        ########## up ##########
        for d in results[4]:
            tdir_up = d["edir"]
            tdif_up = d["edn"]
            data["tdir_up"].append(tdir_up)
            data["tdif_up"].append(tdif_up)

        ########## Spherical albedo ##########
        numerator1 = albedo1 * (np.array(data["rho2"]) - np.array(data["rho0"]))
        numerator2 = albedo2 * (np.array(data["rho1"]) - np.array(data["rho0"]))
        denominator = albedo2 * albedo1 * (np.array(data["rho2"]) - np.array(data["rho1"]))
        spherical_albedo = (numerator1 - numerator2) / denominator
        data["spherical_albedo"] = spherical_albedo

        ########## tdir, tdif, kurudz conversion, edir, edif ##########
        kurudz = np.genfromtxt(libradtran_kurudz_full_path)
        kurudz_filtered = kurudz[(490 <= kurudz[:, 0])
                                 & (kurudz[:, 0] <= 810), :]
        f_kurudz = np.interp(data["wavelength"],
                             kurudz_filtered[:, 0], kurudz_filtered[:, 1])
        edir = np.array(data["tdir_down"]) * f_kurudz * np.cos(sza_radian)
        edif = np.array(data["tdif_down"]) * f_kurudz * np.cos(sza_radian)
        data["edir"] = edir.tolist()
        data["edif"] = edif.tolist()

        ########## Path radiance ##########
        lp0 = np.array(data["rho0"]) * f_kurudz * np.cos(sza_radian) * (1/np.pi)
        data["path_rad"] = lp0.tolist()
        
        # Create a dataframe
        df = pd.DataFrame(data)
        
        # Save df DataFrame into JSON file
        output_file_name = f"{afglms_file_name}_sza{sza}_vza{vza}_phi{phi}_phi0{phi0}_alt{alt}_tau{tau}.json"
        df.to_json(output_file_loc_JSON + output_file_name, orient="columns", indent=4)
    
    save_output_params(results_list)