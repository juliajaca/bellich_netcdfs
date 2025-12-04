# "C:\Users\Julia\Nextcloud\Datos_MM_Art_2025\datasets_ncFormat\Biogeochemical\chlorophyll\EUROGEL\EUROGEL_CHL.nc"

# %%
from netCDF4 import Dataset, stringtochar
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat

# %% COMPROBACION

path = "C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/nutrients/phosphate/EUROGEL/EUROGEL_PO4.nc"  #segundo archivo
path= 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/chlorophyll/IEO_JC/IEO_JC_CHL.nc'
path= 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/chlorophyll/IEO_JM/IEO_JM_CHL.nc'
dataset = Dataset(path, "r")
# %%
print("\n🔹 Dimensiones:")
for var_name in dataset.dimensions:
    print(f"\nDimension: {var_name}")
# %%
vars_nc = {}
print("\n🔹 Atributos de las Variables:")
for var_name in dataset.variables:
    var = dataset.variables[var_name]  # <- accede al objeto variable
    print(f"\nVariable: {var_name}")
    print(f"  Dimensiones: {var.dimensions}")
    print(f"  Shape: {var.shape}")
    print(f" Tipo de dato: {dataset.variables[var_name].dtype}")
    vars_nc[var_name] = dataset.variables[var_name][:]
    for attr in var.ncattrs():
        print(f"  {attr}: {var.getncattr(attr)}")
# %%
print("\n🔹 Atributos Globales:")
for attr in dataset.ncattrs():
    print(f"{attr}: {dataset.getncattr(attr)}")

# %%
for variable, valor in vars_nc.items():
    print(variable, ' tiene ' ,len(valor), ' valores')
    print(valor[0:])
    print('--------------')

# %%