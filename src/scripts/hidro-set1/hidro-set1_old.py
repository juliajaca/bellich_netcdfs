# %%
from netCDF4 import Dataset, stringtochar
import numpy as np
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
from generar_txt import generar_txt

# %%
# nombre_fichero = 'nivel_mar_venom_junio25'

# data0 = pd.read_csv('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/hidrodinamica/set1/Dades-data-2025-04-25 12_39_23.csv', 
#                     sep=',', 
#                     header=0, names=['Fecha', 'imei','bat', 'dis', 'pressio', 'temp', 'write_sd' ],  parse_dates= ['Fecha'] ) 
data0 = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/hidrodinamica/set1/Dades-data-2025-05-15 16_12_06.xlsx',header = None  )
# %%
data1 = data0[0].str.split(",", expand=True)
data1.columns = data1.iloc[0]  # Usa la fila 0 como encabezado
data1.columns = data1.columns.map(str).str.strip(' "\'')
data1 = data1[1:]              # Elimina esa fila del cuerpo de datos

# 4. (Opcional) Resetear índices
data1 = data1.reset_index(drop=True)
data1["Time"] = pd.to_datetime(data1["Time"])
data = data1.sort_values(by=["Time"]).reset_index(drop=True)
data = data.sort_values("Time").drop_duplicates(subset=["Time"]).reset_index(drop=True)
# Convertir la columna a datetime

# Definir la fecha límite
fecha_limite = pd.Timestamp("2024-01-11")

# Filtrar: quedarnos solo con fechas posteriores
data = data[data["Time"] > fecha_limite].copy()
# %%
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (data.Time- epoch) / pd.Timedelta(days=1)
# %%
nombre_fichero= 'VENOM_SSH'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Hydrodynamics/sea_level/VENOM/'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title=f'{nombre_fichero}'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id ='VENOM'
ncfile.project = 'VENOM'; ncfile.source = 'In-situ data collection'; ncfile.Conventions = 'CF-1.8'
ncfile.comments = "Sea level estimated as anomaly from median using a downward-looking sensor: sea_level = -(distance - median(distance))."

# %%
# crear dimensiones
ncfile.createDimension('time', len(dias_desde_1970))
for dim in ncfile.dimensions.items():
    print(dim)

# %%
time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values  # Se asigna directamente
# %%
lat_var = ncfile.createVariable('latitude', np.float64, )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var.grid_mapping = "crs"
lat_var[:] = 37.817740

lon_var = ncfile.createVariable('longitude', np.float64, )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] = -0.783372

# Crear variable con atributo _FillValue al momento de su creación
data["distancia"] = pd.to_numeric(data["distancia"], errors="coerce")
# Calcular anomalía con manejo de NaN
valores_con_nan = -(data["distancia"] - np.nanmedian(data["distancia"]))
# %%
valores_con_nan[np.isnan(valores_con_nan)] = -9999  # Sustituir NaN por el valor faltante
# %%
value_var = ncfile.createVariable('sea_level', np.float64, ('time',), )
value_var.units = 'mm'
value_var.standard_name = 'sea_surface_height_above_mean_sea_level'
value_var.long_name = 'Sea surface height above sea level'
# Convertir columna a numérica (por si hay strings)

# Asignar a la variable NetCDF
value_var[:] = valores_con_nan
value_var.missing_value = -9999  # Opcional, pero útil para software que no reconoce _FillValue
value_var.grid_mapping = "crs"


crs = ncfile.createVariable("crs", "i")  # Dummy variable for coordinate reference system
crs.grid_mapping_name = "latitude_longitude"
crs.projection = "Geodetic"
crs.long_name = "WGS 84 / Geographic coordinates (EPSG:4326)"
crs.epsg_code = "EPSG:4326"
crs.semi_major_axis = 6378137.0
crs.inverse_flattening = 298.257223563
crs.comment = "Geographic coordinates are referenced to WGS 84 (EPSG:4326) in decimal degrees."

ncfile.close()
# %%
#  COPIAR
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Hydrodynamics/sea_level/VENOM/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')


# %% COMPROBACION
dataset = Dataset(f'{path}{nombre_fichero}.nc', "r")
print(dataset.variables.keys())  # Ver las variables en el archivo

print("\n🔹 Atributos de las Variables:")
for var_name in dataset.variables:
    print(f"\nVariable: {var_name}")
    for attr in dataset.variables[var_name].ncattrs():
        print(f"  {attr}: {dataset.variables[var_name].getncattr(attr)}")
print("\n🔹 Atributos Globales:")
for attr in dataset.ncattrs():
    print(f"{attr}: {dataset.getncattr(attr)}")

# %%
# Leer las variables
tiempo = dataset.variables["time"][:]  # Días desde 1970
dist = dataset.variables["sea_level"][:]    #  
print(tiempo); print('-----------------')
print(dist); print('-----------------')

sea_level_raw = dataset.variables["sea_level"][:].data  # Esto ignora la máscara
print("→ Nº de -9999:", np.sum(sea_level_raw == -9999)) #los

# %%
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()
# %% PLOT
fechas = np.array(fechas)  # convertir a numpy array para indexar
orden = np.argsort(fechas)  # obtener índices ordenados

fig, axs = plt.subplots(1, 1, figsize=(15, 10), sharex=True)

axs.plot(fechas, dist, color='green')
axs.set_title('sea_level')
axs.set_ylabel('sea_level')
axs.grid(True)
axs.legend()
plt.tight_layout()
plt.savefig(f'{path}{nombre_fichero}.pdf', format='pdf')
plt.show()

# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')
# %%
"""
 ____  ____  _____ ____  ____  _  ____
/  __\/  __\/  __// ___\/ ___\/ \/  _ \
|  \/||  \/||  \  |    \|    \| || / \|
|  __/|    /|  /_ \___ |\___ || || \_/|
\_/   \_/\_\\____\\____/\____/\_/\____/

"""
nombre_fichero= 'VENOM_ATM_MSLP'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Atmospheric/MLSP/VENOM/'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title=f'{nombre_fichero}'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id ='VENOM'
ncfile.project = 'VENOM'; ncfile.source = 'In-situ data collection'; ncfile.Conventions = 'CF-1.8'


# %%
# crear dimensiones
ncfile.createDimension('time', len(dias_desde_1970))
for dim in ncfile.dimensions.items():
    print(dim)

# %%
time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values  # Se asigna directamente

# %%
lat_var = ncfile.createVariable('latitude', np.float64, )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var.grid_mapping = "crs"
lat_var[:] = 37.817740

lon_var = ncfile.createVariable('longitude', np.float64, )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] = -0.783372

value_var = ncfile.createVariable('air_pressure', np.float64, ('time',))
value_var.units= 'hPa'
value_var.standard_name = 'Air_pressure_at_sea_level'
value_var.long_name= 'Mean sea level atmospheric pressure'
value_var.missing_value = -9999 
value_var.grid_mapping = "crs"

data["pressio"] = pd.to_numeric(data["pressio"], errors="coerce")
valores_con_nan = data.pressio
valores_con_nan[np.isnan(valores_con_nan)] = -9999 
value_var[:] = valores_con_nan

crs = ncfile.createVariable("crs", "i")  # Dummy variable for coordinate reference system
crs.grid_mapping_name = "latitude_longitude"
crs.projection = "Geodetic"
crs.long_name = "WGS 84 / Geographic coordinates (EPSG:4326)"
crs.epsg_code = "EPSG:4326"
crs.semi_major_axis = 6378137.0
crs.inverse_flattening = 298.257223563
crs.comment = "Geographic coordinates are referenced to WGS 84 (EPSG:4326) in decimal degrees."

ncfile.close()
# %% COMPROBACION
dataset = Dataset(f'{path}{nombre_fichero}.nc', "r")
print(dataset.variables.keys())  # Ver las variables en el archivo

print("\n🔹 Atributos de las Variables:")
for var_name in dataset.variables:
    print(f"\nVariable: {var_name}")
    for attr in dataset.variables[var_name].ncattrs():
        print(f"  {attr}: {dataset.variables[var_name].getncattr(attr)}")
print("\n🔹 Atributos Globales:")
for attr in dataset.ncattrs():
    print(f"{attr}: {dataset.getncattr(attr)}")

# %%
# Leer las variables
tiempo = dataset.variables["time"][:]  # Días desde 1970
dist = dataset.variables["air_pressure"][:]    #  
print(tiempo); print('-----------------')
print(dist); print('-----------------')

# %%
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()
# %% PLOT
fechas = np.array(fechas)  # convertir a numpy array para indexar
orden = np.argsort(fechas)  # obtener índices ordenados

fig, axs = plt.subplots(1, 1, figsize=(15, 10), sharex=True)

axs.plot(fechas, dist, color='purple', )
axs.set_title('air pressure')
axs.set_ylabel('air pressure')
axs.grid(True)
axs.legend()
plt.tight_layout()
plt.savefig(f'{path}{nombre_fichero}.pdf', format='pdf')
plt.show()

# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Atmospheric/MLSP/VENOM/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')
# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')

# %% SE DECIDE QUITAR PORQUE NO ES FIABLE 
# """
#   _______ .____  __   __ .___  .____  .___          .
#  '   /    /      |    |  /   \ /      /   \   ___  _/_   ,   . .___    ___
#      |    |__.   |\  /|  |,_-' |__.   |__-'  /   `  |    |   | /   \  /   `
#      |    |      | \/ |  |     |      |  \  |    |  |    |   | |   ' |    |
#      /    /----/ /    /  /     /----/ /   \ `.__/|  \__/ `._/| /     `.__/|

# """
# nombre_fichero= 'VENOM_TEMP_AIR'
# path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Atmospheric/temperature/VENOM/'
# ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
# print(ncfile)

# # %% CREAR ATRIBUTOS GLOBALES
# ncfile.title=f'{nombre_fichero}'
# ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
# ncfile.domain= 'Mar menor coastal lagoon, Spain'
# ncfile.dataset_id ='VENOM'
# ncfile.project = 'VENOM'; ncfile.source = 'In-situ data collection'; ncfile.Conventions = 'CF-1.8'

# # %%
# # crear dimensiones
# ncfile.createDimension('time', len(dias_desde_1970))

# for dim in ncfile.dimensions.items():
#     print(dim)

# # %%
# time_var = ncfile.createVariable('time', np.float64, ('time',))
# time_var.units = "days since 1970-01-01 00:00:0"
# time_var.standard_name = "time"
# time_var.calendar = 'gregorian'
# time_var[:] = dias_desde_1970.values  # Se asigna directamente

# # %%
# lat_var = ncfile.createVariable('latitude', np.float64, )
# lat_var.units = 'degrees_north'
# lat_var.standard_name = "latitude"
# lat_var.grid_mapping = "crs"
# lat_var[:] = 37.817740

# lon_var = ncfile.createVariable('longitude', np.float64, )
# lon_var.units = 'degrees_east'
# lon_var.standard_name = "longitude"
# lon_var.grid_mapping = "crs"
# lon_var[:] = -0.783372

# value_var = ncfile.createVariable('air_temperature', np.float64, ('time',))
# value_var.standard_name = 'Air_temperature'
# value_var.long_name = 'Mean air temperature'
# value_var.units= 'degree_Celsius'
# data["temperatura"] = pd.to_numeric(data["temperatura"], errors="coerce")
# value_var[:] = data.temperatura

# value_var = ncfile.createVariable('air_pressure', np.float64, ('time',))
# value_var.standard_name = 'Air_pressure_at_sea_level'
# value_var.units= 'hPa'
# value_var.grid_mapping = "crs"
# value_var.missing_value = -9999 
# data["pressio"] = pd.to_numeric(data["pressio"], errors="coerce")
# valores_con_nan = data.pressio
# valores_con_nan[np.isnan(valores_con_nan)] = -9999 
# value_var[:] = valores_con_nan


# crs = ncfile.createVariable("crs", "i")  # Dummy variable for coordinate reference system
# crs.grid_mapping_name = "latitude_longitude"
# crs.projection = "Geodetic"
# crs.long_name = "WGS 84 / Geographic coordinates (EPSG:4326)"
# crs.epsg_code = "EPSG:4326"
# crs.semi_major_axis = 6378137.0
# crs.inverse_flattening = 298.257223563
# crs.comment = "Geographic coordinates are referenced to WGS 84 (EPSG:4326) in decimal degrees."

# ncfile.close()
# # %% COMPROBACION
# dataset = Dataset(f'{path}{nombre_fichero}.nc', "r")
# print(dataset.variables.keys())  # Ver las variables en el archivo

# print("\n🔹 Atributos de las Variables:")
# for var_name in dataset.variables:
#     print(f"\nVariable: {var_name}")
#     for attr in dataset.variables[var_name].ncattrs():
#         print(f"  {attr}: {dataset.variables[var_name].getncattr(attr)}")
# print("\n🔹 Atributos Globales:")
# for attr in dataset.ncattrs():
#     print(f"{attr}: {dataset.getncattr(attr)}")

# # %%
# # Leer las variables
# tiempo = dataset.variables["time"][:]  # Días desde 1970
# dist = dataset.variables["air_temperature"][:]    #  
# print(tiempo); print('-----------------')
# print(dist); print('-----------------')

# # %%
# fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
# dataset.close()
# # %% PLOT
# fechas = np.array(fechas)  # convertir a numpy array para indexar
# orden = np.argsort(fechas)  # obtener índices ordenados

# fig, axs = plt.subplots(1, 1, figsize=(15, 10), sharex=True)

# axs.plot(fechas, dist, color='red',)
# axs.set_title('air temp')
# axs.set_ylabel('air temps')
# axs.grid(True)
# axs.legend()
# plt.tight_layout()
# plt.savefig(f'{path}{nombre_fichero}.pdf', format='pdf')
# plt.show()

# # %%
# generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')
# #  
# # %%
# ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Atmospheric/temperatura/VENOM/'
# shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')