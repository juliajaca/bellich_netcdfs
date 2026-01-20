# %%
from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
from generar_txt import generar_txt
# %%
data  = pd.read_excel('C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/BioGeoquimicas/IEO_SAT_GomezJakobsen/salida3.xlsx',sheet_name="CHLORO")

data = data.sort_values(by=['anio', 'mes','dia']).reset_index(drop=True)
print(data.head())
data["fecha"] = pd.to_datetime(
    dict(year=data["anio"], month=data["mes"], day=data["dia"])
)
data = data.dropna(subset=data.columns[3:29]
, how="all")
# %%
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (data.fecha - epoch) / pd.Timedelta(days=1)
# %%
nombre_fichero = 'BELICH_SAT_BELA_TS_CHL'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/chlorophyll/BELICH_SAT_BELA_TS/'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)


# %% CREAR ATRIBUTOS GLOBALES
ncfile.title= nombre_fichero
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'BELICH_SAT_BELA_TS'
ncfile.project = 'BELICH'; ncfile.source = 'Satellite data'; ncfile.Conventions = 'CF-1.8'
ncfile.comment = 'Corrected by BELA algorithm'

# %%
# crear dimensiones
max_length_param = len("C10")
estaciones = data.columns[3:29].astype(str).tolist()
# longitud máxima automáticamente
max_length = max(len(s) for s in estaciones)
# array numpy compatible con netCDF
estaciones_np = np.array(estaciones, dtype=f"S{max_length}")

# Asiganar posiciones a las estaciones

df_pos = pd.read_excel('C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/BioGeoquimicas/IEO_SAT_GomezJakobsen/salida3.xlsx',sheet_name="POSI")
df_pos_filtrado = df_pos[df_pos["nom"].isin(estaciones)]
df_pos_ordenado = (
    df_pos_filtrado
    .set_index("nom")
    .loc[estaciones]
    .reset_index()
)
latitudes = df_pos_ordenado["lat"].to_numpy()
longitudes = df_pos_ordenado["lon"].to_numpy()

print(latitudes.shape, longitudes.shape)
# %%

ncfile.createDimension('time', len(dias_desde_1970))
ncfile.createDimension('unit_char_len', max_length_param)
ncfile.createDimension('station', len(estaciones))
for dim in ncfile.dimensions.items():
    print(dim)

# %%
time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values  # Se asigna directamente

lat_var = ncfile.createVariable('latitude', np.float64,('station') )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var.grid_mapping = "crs"
lat_var[:] = latitudes

lon_var = ncfile.createVariable('longitude', np.float64,('station') )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] = longitudes

# %%
cols_est = data.columns[3:29]# columnas de estaciones
# array 2D (fecha, estacion)
df_chl = data[cols_est].to_numpy()
value_var = ncfile.createVariable('chlorophyll', np.float64, ('time','station'))
value_var.units= 'mg m-3'
value_var.standard_name = "mass_concentration_of_chlorophyll_a_in_sea_water"
value_var.long_name = 'Chlorophyll-a Concentration in Sea Water'
value_var.cell_methods = "time: mean"
value_var.missing_value = -9999
value_var.grid_mapping = "crs"
value_var[:,:] = df_chl


parameter_var = ncfile.createVariable('station_name', 'S1', ('station', 'unit_char_len'))
parameter_var.long_name = 'station'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)

crs = ncfile.createVariable("crs", "i")  # Dummy variable for coordinate reference system
crs.grid_mapping_name = "latitude_longitude"
crs.projection = "Geodetic"
crs.long_name = "WGS 84 / Geographic coordinates (EPSG:4326)"
crs.epsg_code = "EPSG:4326"
crs.semi_major_axis = 6378137.0
crs.inverse_flattening = 298.257223563
crs.comment = "Geographic coordinates are referenced to WGS 84 (EPSG:4326) in decimal degrees."
# %%
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
unit = dataset.variables["station_name"][:]    #  
value = dataset.variables["chlorophyll"][:] 
print(tiempo); print('-----------------')
print(unit); print('-----------------')
print(value); print('-----------------')
# %%
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()
# %% PLOT
fechas = np.array(fechas)  # convertir a numpy array para indexar
orden = np.argsort(fechas)  # obtener índices ordenados

# Reordenar fechas y nitrato
fechas_ordenadas = fechas[orden]
nitrato_ordenado = value[orden, :]

fig, axes = plt.subplots( nrows=13, ncols=2, figsize=(20, 40), sharex=True)
axes = axes.flatten()

for i in range(26):
    axes[i].plot(fechas_ordenadas, nitrato_ordenado[:, i], marker='o', label=estaciones[i])
    axes[i].grid(True)
    axes[i].legend()

    if i % 2 == 0:
        axes[i].set_ylabel("clorofila")
    else:
        axes[i].set_ylabel("")

axes[-1].set_xlabel('Fecha')
plt.tight_layout()
plt.show()

# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}.txt')
# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Biogeochemical/chlorophyll/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')
# %%