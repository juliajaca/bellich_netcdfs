# %%
from netCDF4 import Dataset, stringtochar
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat

# %%
# Cargar el archivo .mat
data = loadmat('C:/Users/Julia/Documents/VSCODE_BELLICH/src/scripts/_REVISION_LILIA/salinidad/set1/CTD profiles campañas 2016-17/CTD profiles campañas 2016-17/MM201611_ctd_profiles.mat')

# Mostrar las variables disponibles en el archivo
print(data.keys())


# %% COMPROBACION
dataset = Dataset(f'C:/Users/Julia/Documents/VSCODE_BELLICH/src/scripts/_REVISION_LILIA/salinidad/set1/set1_feb2017_IEO/set1_feb2017_SAL_IEO.nc', "r")
print(dataset.variables.keys())  # Ver las variables en el archivo

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
    print(valor[0:10])
    print('--------------')
# %%

fechas = pd.to_datetime(vars_nc['time'], origin="1970-01-01", unit="D")
fechas = np.array(fechas) 

time = dataset.variables['time'][:]
temp = dataset.variables['seawater_temperature'][:]
# %%
# %%
dataset.close()

# %%
idx_tiempo = np.where(fechas == np.datetime64('2004-12-01'))[0]

# Encuentra el índice de profundidad más cercano a -0.5 m
idx_depth = np.where(vars_nc['depth'].data == 0)[0]

# %%
temp_slice = vars_nc['seawater_temperature'][idx_tiempo, idx_depth, :, :]  # shape: (lat, lon)
temp_slice = np.squeeze(temp_slice)
# Crea el grid de coordenadas
lon_grid, lat_grid = np.meshgrid(vars_nc['longitude'], vars_nc['latitude'])

# Aplana todo para scatterplot
lon_flat = lon_grid.flatten()
lat_flat = lat_grid.flatten()
temp_flat = temp_slice.flatten()

# Haz el scatterplot
plt.figure(figsize=(8, 6))
sc = plt.scatter(lon_flat, lat_flat, c=temp_flat, cmap='viridis', s=30)
plt.colorbar(sc, label='Temperatura del agua (°C)')
plt.title('Temperatura a -0 m el 2004-12-01')
plt.xlabel('Longitud')
plt.ylabel('Latitud')
plt.grid(True)
plt.show()
# %%
# %% PERFIL DE TEMPERATURA
idx_tiempo = np.where(fechas == np.datetime64('2003-01-22'))[0]
idx_lat = np.where(vars_nc['latitude'].data == [37.7833])
idx_lon = np.where(vars_nc['longitude'].data == [-0.78333])
temp_profile = vars_nc['seawater_temperature'][idx_tiempo, :, idx_lat, idx_lon]
temp_profile = np.squeeze(temp_profile)

# Paso 6: plot del perfil
plt.figure(figsize=(6, 8))
plt.plot(temp_profile[0], vars_nc['depth'], marker='o')
plt.gca().invert_yaxis()  # Invertir eje Y para que la profundidad aumente hacia abajo
plt.xlabel('Temperatura (°C)')
plt.ylabel('Profundidad (m)')
plt.title('Perfil de Temperatura\n8 de marzo de 2022, lat 37.6771, lon -0.7507')
plt.grid(True)
plt.show()

# %% SERIE TEMPORAL
idx_depth = np.where(vars_nc['depth'].data == 1)[0]
temp_timeseries = vars_nc['seawater_temperature'][:, idx_depth, idx_lat, idx_lon]
temp_timeseries = np.squeeze(temp_timeseries)

plt.figure(figsize=(10, 5))
plt.plot(fechas, temp_timeseries, marker='o')
plt.title(f'Serie temporal de temperatura\nLat: lat 37.6771, lon -0.7507 a 1.5 m')
plt.xlabel('Fecha')
plt.ylabel('Temperatura (°C)')
plt.grid(True)
plt.tight_layout()
plt.show()