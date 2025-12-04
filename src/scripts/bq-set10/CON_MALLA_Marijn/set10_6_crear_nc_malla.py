# %%
from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import pandas as pd
import subprocess
import shutil
import os
from datetime import datetime
import matplotlib.pyplot as plt
import sys
sys.path.append("../../")
from generar_txt import generar_txt


# %%
data = pd.read_csv('C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/BioGeoquimicas/2023-2025_BELA_clorofila_de_marijn_con_malla_CSV/20230105_BELA.csv', decimal='.', sep=',')

folder = "C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/BioGeoquimicas/2023-2025_BELA_clorofila_de_marijn_con_malla_CSV/"

# Listamos todos los archivos csv
files = sorted([f for f in os.listdir(folder) if f.endswith(".csv")])

dates_list = []
data_list = []

# %%
for f in files:
    # Extraemos la fecha del nombre (primeros 8 dígitos)
    date_str = f[:8]
    date = datetime.strptime(date_str, "%Y%m%d")
    print(date)
    print(type(date))

    dates_list.append(date)

    # Leemos el csv
    df = pd.read_csv(os.path.join(folder, f))

    # Guardamos los valores de lectura (la "capa" para esa fecha)
    data_list.append(df["chl"].values)

# %%
# Convertimos fechas a array numpy
dates = np.sort(pd.Series(dates_list).unique())
fechas_unicas = np.array(dates)
fechas_unicas_ts = pd.to_datetime(fechas_unicas)
# Calcular la diferencia en días
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (fechas_unicas_ts - epoch) / pd.Timedelta(days=1)

# Obtenemos lat/lon del primer archivo (suponiendo todos iguales)
sample = pd.read_csv(os.path.join(folder, files[0]))
lats = sample["latitude"].values
lons = sample["longitude"].values

# Convertimos lista de arrays en un array 3D: time × points
data_2d = np.array(data_list)  

# %%

path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/chlorophyll/BELICH_SAT_BELA_MAPS/'
nombre_fichero= 'BELICH_SAT_BELA_MAPS'
# %%
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF4_CLASSIC')
# NETCDF3_CLASSIC siempre prealoca en disco toda la variable al crearla.
# Si tus dimensiones son grandes, esa preasignación puede tardar minutos.
# asi que creamos un netcdf4
print(ncfile)

# %%
ncfile.title=nombre_fichero
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'BELICH_SAT_BELA_MAPS'
ncfile.project = 'Belich'
ncfile.source = 'Copernicus'
ncfile.Conventions = 'CF-1.8'

# Crear dimensiones
ncfile.createDimension("time", len(fechas_unicas)) 
ncfile.createDimension("latitude", len(np.unique(lats))) 
ncfile.createDimension("longitude", len(np.unique(lons)))

for dim in ncfile.dimensions.items():
    print(dim)
# %%
time_var = ncfile.createVariable('time', np.float64, ('time'))
time_var.units= "days since 1970-01-01 00:00:0"
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
time_var[:] = dias_desde_1970.values

lat_var = ncfile.createVariable('latitude', np.float64, ('latitude',))
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var.grid_mapping = "crs"
lat_var[:] = np.unique(lats)

lon_var = ncfile.createVariable('longitude', np.float64, ('longitude',))
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] = np.unique(lons)

# %%
chl_var = ncfile.createVariable('chlorophyll', np.float64, ("time", "latitude", "longitude"), fill_value=None)
chl_var.units = 'ug L-1'
chl_var.standard_name = 'mass_concentration_of_chlorophyll_a_in_sea_water'
chl_var.long_name = 'Surface chlorophyll concentration from satellite'
chl_var.missing_value = -9999
chl_var.grid_mapping = "crs"
chl_var.comment = 'BELA chlorophyll concentration algorithm used to estimate CHL maps.'

# %%
# Escribir datos en el NetCDF (fila por fila)
unique_lats = np.unique(lats)
unique_lons = np.unique(lons)
lat_to_idx = {lat: i for i, lat in enumerate(unique_lats)}
lon_to_idx = {lon: i for i, lon in enumerate(unique_lons)}

# Mapa 2D: par lat-lon → celda de la matriz
coord_to_index = {
    (lat, lon): (lat_to_idx[lat], lon_to_idx[lon])
    for lat, lon in zip(lats, lons)
}

for i, fecha in enumerate(fechas_unicas_ts):
    print(fecha)
    sub = data_2d[i]   # DataFrame chl de esa fecha
    chl_matrix = np.full((len(unique_lats), len(unique_lons)), -9999,dtype= np.float64)

    for lat, lon, chl in zip(lats, lons, sub):
        # índice rápido O(1)
        i_lat, i_lon = coord_to_index[(lat, lon)]
        # print(i_lat, i_lon)
        chl_matrix[i_lat, i_lon] = chl

    chl_var[i, :, :] = chl_matrix

# %%
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
lat = dataset.variables["latitude"][:] 
lon = dataset.variables["longitude"][:] 
value = dataset.variables["chlorophyll"][:,:,:].data
print(tiempo); print('-----------------')
print(lat); print('-----------------')
print(lon); print('-----------------')
print(value); print('-----------------')
valores_validos = value[0, :, :][value[0, :, :] != -9999]
print(valores_validos)


print(valores_validos)

# %%
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()

# %% PLOT
#  --- Seleccionar índices para 5 días ---
nt = len(tiempo)
indices_dias = [0, nt//4, nt//2, 3*nt//4, nt-1]

for idx in indices_dias:
    plt.figure(figsize=(8,6))
    
    # Datos del día
    data_dia = value[idx,:,:].astype(np.float64)
    
    lon2d, lat2d = np.meshgrid(lon, lat)  # shapes: (n_lat, n_lon)

    # Reemplazar valores missing
    data_dia[data_dia == -9999] = np.nan
    
    # Extraer coordenadas válidas
    valid = ~np.isnan(data_dia)
    lat_valid = lat2d[valid]
    lon_valid = lon2d[valid]
    chl_valid = data_dia[valid]
    
    # Scatter plot
    plt.scatter(lon_valid, lat_valid, c=chl_valid, cmap='viridis', s=20)  # <--- s controla el tamaño
    plt.colorbar(label='Chlorophyll [ug/L]')
    
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title(f'Chlorophyll - Día {idx} (time={tiempo[idx]})')
    
    plt.tight_layout()
    plt.show()


# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}.txt')
# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Biogeochemical/chlorophyll/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')
# %%
