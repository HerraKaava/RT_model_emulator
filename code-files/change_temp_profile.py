# In[]:
# ########################################################################
# Packages
# ########################################################################

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

# In[]:
# ########################################################################
# Reading from a text file with NumPy
# ########################################################################

# Specify file path that you want to read in
file_path = "/home/jaminkiukkonen/libRadtran-2.0.5/data/atmmod/afglms.dat"

# Read in the file as a NumPy array
arr = np.genfromtxt(file_path, dtype='float')

# The third column contains the temperature profile
print(arr[:, 2].shape)    # (50,)

# In[]:
# ########################################################################
# Read in the seasonal temperature profiles
# ########################################################################

def read_netCDF4_data(filename):
    with Dataset(filename) as rootgrp:
        t = rootgrp["data"][:, :]
        return t

# Seasonal temperature profiles
Q1_profiles = read_netCDF4_data("/home/jaminkiukkonen/Desktop/work/SummerProject/data/Q1_profiles.nc")
Q2_profiles = read_netCDF4_data("/home/jaminkiukkonen/Desktop/work/SummerProject/data/Q2_profiles.nc")
Q3_profiles = read_netCDF4_data("/home/jaminkiukkonen/Desktop/work/SummerProject/data/Q3_profiles.nc")
Q4_profiles = read_netCDF4_data("/home/jaminkiukkonen/Desktop/work/SummerProject/data/Q4_profiles.nc")

print(Q1_profiles.shape)    # (18, 19)
print(Q2_profiles.shape)    # (18, 19)
print(Q3_profiles.shape)    # (18, 19)
print(Q4_profiles.shape)    # (18, 19)

# In[]:
# ########################################################################
# Interpolate the temperature profiles to match the length of the
# midlatitude summer temperature profile
# ########################################################################

new_pressure_levels = arr[:, 1]    # The 2nd col contains the pressure levels
original_pressure_levels = np.array([1, 5, 10, 30, 50, 70, 100, 125, 175, 225, 300, 400, 500, 600, 700, 800, 900, 950, 1000])
temperatures = Q1_profiles[0, :]
interpolated_temperatures = np.interp(new_pressure_levels, original_pressure_levels, temperatures)

# Visualize
plt.plot(interpolated_temperatures, new_pressure_levels, c="red")
plt.plot(temperatures, original_pressure_levels, c="black", linestyle="--")
plt.gca().invert_yaxis()

# In[]:
# ########################################################################
# Changing the temperature profile of the arr (midlatitude summer file)
# ########################################################################

arr[:, 1] = interpolated_temperatures

# In[]:
# ########################################################################
# Writing a text file from a NumPy array

col_names = " z(km), p(mb), T(K), air(cm-3), o3(cm-3), o2(cm-3), h2o(cm-3), co2(cm-3), no2(cm-3)"
np.savetxt("test_arr.txt", arr, header=col_names, delimiter=",", fmt='%20.6e', comments="#")