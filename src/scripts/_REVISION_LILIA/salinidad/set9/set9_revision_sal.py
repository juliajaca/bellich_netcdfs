# %%
from netCDF4 import Dataset, stringtochar
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat

# %% COMPROBACION
dataset = Dataset(f'C:/Users/Julia/Documents/VSCODE_BELLICH/src/scripts/_REVISION_LILIA/salinidad/set9/set9_SAL_IEO.nc', "r")

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

# %%
plt.figure(figsize=(10, 5))
plt.plot(vars_nc['time'], vars_nc['seawater_salinity'][:,1], marker='o')
plt.xlabel('Fecha')
plt.ylabel('salinidad')
plt.grid(True)
plt.tight_layout()
plt.show()
# %%
plt.figure(figsize=(10, 5))
plt.plot(vars_nc['time'], vars_nc['seawater_salinity'][:,0], marker='o')
plt.xlabel('Fecha')
plt.ylabel('salinidad')
plt.grid(True)
plt.tight_layout()
plt.show()

# %%
plt.figure(figsize=(10, 5))
plt.plot(vars_nc['time'], vars_nc['seawater_salinity'][:,2], marker='o')
plt.xlabel('Fecha')
plt.ylabel('salinidad')
plt.grid(True)
plt.tight_layout()
plt.show()
# %%
dataset.close()