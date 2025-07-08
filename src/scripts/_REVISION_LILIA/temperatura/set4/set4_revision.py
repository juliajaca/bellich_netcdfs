# %%
from netCDF4 import Dataset, stringtochar
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# %% COMPROBACION
dataset = Dataset(f'C:/Users/Julia/Documents/VSCODE_BELLICH/src/scripts/_REVISION_LILIA/temperatura/set4/set4_IEO_DMMEM/set4_TEMP_IEO_DMMEM.nc', "r")
print(dataset.variables.keys())  # Ver las variables en el archivo

print("\n🔹 Dimensiones:")
for var_name in dataset.dimensions:
    print(f"\nDimension: {var_name}")

 
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

print("\n🔹 Atributos Globales:")
for attr in dataset.ncattrs():
    print(f"{attr}: {dataset.getncattr(attr)}")

# %%

fechas = pd.to_datetime(vars_nc['time'], origin="1970-01-01", unit="D")
fechas = np.array(fechas) 
dataset.close()



# %%


