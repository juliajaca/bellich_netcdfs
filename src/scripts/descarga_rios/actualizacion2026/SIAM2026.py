# %%
from netCDF4 import Dataset, stringtochar
import numpy as np
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("../../")
from generar_txt import generar_txt
# %%
n_dataset ='SIAM'
nombre_fichero = 'SIAM_PR'
data = pd.read_csv('C:/Users/Julia/Documents/VSCODE_BELLICH/src/scripts/descarga_rios/actualizacion2026/datos/SIAM/SIAM_rainfall.csv',  sep=',', parse_dates= ['DateTime'],) #son todo strings
data.groupby(['ID_station']).size()

data['DateTime'] = data['DateTime'].astype(str)  # asegurar que sean strings
data['DateTime'] = data['DateTime'].apply(lambda x: x if ':' in x else x + ' 00:00:00') #hay agunas estaciones que a las 0:00 no pone la hora
# %%

data['DateTime'] = pd.to_datetime(data['DateTime'], errors='coerce',)
# %% hay otras estaciones que las 0:00 pone el dia anterior 23:50
# Crear máscara de filas con hora 23:50
mask_2350 = data['DateTime'].dt.hour.eq(23) & data['DateTime'].dt.minute.eq(50)

# Pasarlas al día siguiente 00:00
data.loc[mask_2350, 'DateTime'] = data.loc[mask_2350, 'DateTime'] + pd.Timedelta(minutes=10)

data = data.sort_values(by=["ID_station", "DateTime"]).reset_index(drop=True)

unique_times = data['DateTime'].drop_duplicates().sort_values().reset_index(drop=True)
epoch = pd.Timestamp('1970-01-01')
time_values = (unique_times - epoch) / pd.Timedelta(days=1)

# %%
estaciones = data['ID_station'].unique()
estaciones_np = np.array(estaciones, dtype=f'S{max(len(s) for s in estaciones)}')
max_length_param = max(len(s) for s in estaciones)

coordenadas = pd.read_csv('C:/Users/Julia/Documents/VSCODE_BELLICH/src/scripts/descarga_rios/actualizacion2026/datos/SIAM/SIAM_stations.csv')
nombres_completos = [coordenadas.iloc[i].Name for i in range(len(coordenadas))]
max_len_nombre = max(len(s) for s in nombres_completos)

muni_completos = [coordenadas.iloc[i].Municipality for i in range(len(coordenadas))]
max_len_muni = max(len(s) for s in muni_completos)

inst_completas = [coordenadas.iloc[i].Date_setup for i in range(len(coordenadas))]
max_len_fecha = max(len(str(f)) for f in inst_completas)  # longitud máxima

# %%
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Atmospheric/precipitation/SIAM/'
path_copia = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Atmospheric/precipitation/'

ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')

ncfile.title = nombre_fichero
ncfile.institution="Servicio de Información Agrometeorológica de Murcia (SIAM; Murcia Agrarian Information System)"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = n_dataset
ncfile.project = n_dataset; ncfile.source = 'In situ data collection'; ncfile.Conventions = 'CF-1.8'
ncfile.comment ='Meteorological stations comply with Standard UNE 176.101:2010'

'''
SIAM Servicio agrometeorológico de Murcia en la cuenca vertiente del Mar Menor, 12 stations, hourly sampling. (dataset: Precip_SIAM_16_24_Mma.xls, Precip_SIAM_16_24_Mmb.xls,  Estac_coord.xls, responsible: Gonzalo González Barberá, period covered: 1/01/2016 to 30/11/2024 pero varia con la estacion, path: Datos_MM_Art_2025\Atmosfericas\GGB, variables measured: hourly precipitation).
'''

all_dates = pd.date_range(start=data['DateTime'].min(),
                          end=data['DateTime'].max(),
                          freq='H')

# Crear pivot table de precipitación
pivot = data.pivot_table(index='DateTime', columns='ID_station', values='Rainfall')

# Reindex para incluir TODAS las horas del rango, poniendo NaN donde no hay datos
pivot = pivot.reindex(all_dates)

# Convertir a numpy array y reemplazar NaN por -9999
valor_array = pivot.to_numpy()
valores_con_nan = valor_array.copy()
valores_con_nan[np.isnan(valores_con_nan)] = -9999

# Variable time
epoch = pd.Timestamp('1970-01-01')
time_values = (all_dates - epoch) / pd.Timedelta(days=1)

ncfile.createDimension('time', len(all_dates))
ncfile.createDimension('station', len(estaciones))
ncfile.createDimension('station_full_strlen', max_length_param)
ncfile.createDimension('name_full_strlen', max_len_nombre)
ncfile.createDimension('municipality_full_strlen', max_len_muni)
ncfile.createDimension('installation_date_strlen', max_len_fecha)

for dim in ncfile.dimensions.items():
    print(dim)
# %%
time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:00"
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
dias_desde_1970 = pd.Series(pivot.index)
dias_desde_1970 = (dias_desde_1970 - epoch) / pd.Timedelta(days=1)
time_var[:] = time_values.values # Se asigna directamente

# Latitud y longitud
lat_var = ncfile.createVariable('latitude', np.float64, ('station',))
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var.grid_mapping = "crs"
lat_var[:] =  [coordenadas.iloc[i].Latitude for i in range(len(coordenadas))]

lon_var = ncfile.createVariable('longitude', np.float64, ('station',))
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] =  [coordenadas.iloc[i].Longitude for i in range(len(coordenadas))]

valores_con_nan = valor_array.copy()
valores_con_nan[np.isnan(valores_con_nan)] = -9999
value_var = ncfile.createVariable('precipitation', np.float64, ('time', 'station'))
value_var.units= 'mm'
value_var.standard_name = 'lwe_thickness_of_precipitation_amount'
value_var.long_name = 'hourly accumulated precipitation'
value_var.cell_methods= 'time: sum'
value_var.missing_value = -9999
value_var.grid_mapping = "crs"
value_var.time_resolution = '1 hour'
value_var[:,:] = valores_con_nan

parameter_var = ncfile.createVariable('station_code', 'S1', ('station', 'station_full_strlen'))
parameter_var.long_name = 'station code'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)

nombre_var = ncfile.createVariable('station_name', 'S1', ('station', 'name_full_strlen'))
nombre_var.long_name = 'station name'
nombre_var._Encoding = 'ascii'
nombre_var[:, :] = stringtochar(np.array(nombres_completos, dtype=f'S{max_len_nombre}'))

muni_var = ncfile.createVariable('station_municipality', 'S1', ('station', 'municipality_full_strlen'))
muni_var.long_name = 'station municipality'
muni_var._Encoding = 'ascii'
muni_var[:, :] = stringtochar(np.array(muni_completos, dtype=f'S{max_len_muni}'))

instal_var = ncfile.createVariable('installation_date', 'S1', ('station', 'installation_date_strlen'))
instal_var.long_name = "station installation date"
instal_var._Encoding = 'ascii'
instal_var[:, :] = stringtochar(np.array(inst_completas, dtype=f'S{max_len_fecha}'))
instal_var.comment = 'Since the time the stations were setup there may be that the station were not recording data because diverse failures. Data referred to these periods are not labelled as missing, but simply they are not in the dataset'

alt = ncfile.createVariable('altitude', np.float64, ('station',))
alt.long_name = "station altitude"
alt.units = 'm'
alt[:] = [coordenadas.iloc[i].Altitude for i in range(len(coordenadas))]

catch = ncfile.createVariable('catchment', np.int8, ('station'))
catch.long_name = "station in Mar Menor Catchment"
catch.flag_values = np.array([0,1], dtype=np.int8)
catch.flag_meanings = "No Yes"
catchment_list = [coordenadas.iloc[i].MarMenor_catchmmet for i in range(len(coordenadas))]
catch[:] = [1 if x=="Y" else 0 for x in catchment_list]
catch.comment = 'the station is within the boundaries of the terrestrial catchment of the Mar Menor (1) or not, but nearby (0)'

crs = ncfile.createVariable("crs", "i") 
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

tiempo = dataset.variables["time"][:]  # Días desde 1970
value = dataset.variables['precipitation'][:] 
station_name = dataset.variables['station_name'][:]
station_code  = dataset.variables['station_code'][:]
station_muni  = dataset.variables['station_municipality'][:]
station_date = dataset.variables['installation_date'][:]
station_altitude = dataset.variables['altitude'][:]
station_c = dataset.variables['catchment'][:]
lat = dataset.variables['latitude'][:]
lon = dataset.variables['longitude'][:]

print(tiempo); print('-----------------')
print(station_code); print('-----------------')
print(value[0:10]); print('-----------------')
print(station_name); print('-----------------')

fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()

# %% PLOT
fechas = np.array(fechas)  # convertir a numpy array para indexar
orden = np.argsort(fechas)  # obtener índices ordenados

# Reordenar fechas y nitrato
fechas_ordenadas = fechas[orden]

fig, axs = plt.subplots(9, 2, figsize=(20,40), sharex=True, sharey=True)
axes = axs.flatten()

for i in range(17):
    axes[i].plot(fechas, value[:, i], marker='o', linestyle='None',markersize=4 )
    axes[i].grid(True)

    axes[i].text(
        0.02, 0.95,
        f"{station_name[i]} - {station_code[i]} - Altitud: {station_altitude[i]} m\n -Catchment: {station_c[i]}- Municipio {station_muni[i]} - instalación {station_date [i]} lat {lat [i]} y lon  {lon [i]}",
        transform=axes[i].transAxes,
        verticalalignment='top',
        horizontalalignment='left',
        fontsize=10,
        bbox=dict(facecolor='white', alpha=0.7, boxstyle="round,pad=0.3")
)    
    if i % 2 == 0:
        axes[i].set_ylabel("precipitacion")
    else:
        axes[i].set_ylabel("")

axes[-1].set_xlabel('Fecha')
plt.tight_layout()
plt.savefig(
    f"{path}SIAM_precipitation_all_stations.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()

# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')
# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Atmospheric/precipitation/SIAM/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')
# %% pruebas
est1 = data.loc[data['ID_station'] == 'AL52']
est1['Rainfall'].plot(style='o', title='Precipitación estación AL52')
# %%
est1 = data.loc[data['ID_station'] == 'CA12']
est1['Rainfall'].plot(style='o', title='Precipitación estación ca12')

# %%
est1 = data.loc[data['ID_station'] == 'TP52']
est1['Rainfall'].plot(style='o', title='Precipitación estación tp52')

# %%
