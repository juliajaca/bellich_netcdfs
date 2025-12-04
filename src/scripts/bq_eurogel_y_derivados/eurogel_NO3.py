# %%
from netCDF4 import Dataset, stringtochar
import numpy as np
import pandas as pd
import shutil
import os
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
from generar_txt import generar_txt
# %%
n_dataset ='EUROGEL'
nombre_fichero = 'EUROGEL_NO3'

data0 = pd.read_excel('C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Physico-chemical/IEO_EUROGEL/060325_resumen_historicos.xlsx',  dtype={  
"valor": "float64", }, parse_dates= ['fecha'],  usecols= ['var', 'estacion', 'fecha','profundidad','valor' , 'muestreo', 'latitud', 'longitud'])

# %%
oxy = data0.loc[(data0['var'] == 'nitrato') & (data0['muestreo'] == 'EUROGEL')]
print(data0.head())
print(f'tiene una longitud de {len(data0)} filas')
oxy = oxy.replace(np.nan, -9999) #reemplazo los nan por -9999 
print(oxy.tail())
# %%
max_length_param = len("RAMBLA")
estaciones = ['E01', 'E02', 'E03', 'E04', 'E05', 'E06', 'E07', 'E08', 'E09',
       'E10', 'E11', 'E12', 'E13', 'E14', 'E15', 'E16', 'E17', 'E18',
       'E19','E20', 'E21', 'E22', 'E23', 'E24', 'RAMBLA' ]
estaciones_np = np.array(estaciones, dtype=f'S{max(len(s) for s in estaciones)}')
depths = oxy.profundidad.unique()

# Creamos una máscara booleana: True si el número es entero
mask_enteros = (depths == np.floor(depths))
# Aplicamos la máscara
depths_enteros = [0,1]
depths_bottom = np.sort(depths[~mask_enteros])

# %%
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/nutrients/nitrate/EUROGEL/'
path_copia = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Biogeochemical/nutrients/nutrients/nitrate/'

data = oxy.sort_values(by=["fecha", "estacion", 'profundidad']).reset_index(drop=True)
print('la longitud antes es', len(data))
data = data.drop_duplicates()
print('la longitud despues es', len(data))
data['fecha'] = pd.to_datetime(data['fecha'], errors='coerce')
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (data.fecha.drop_duplicates() - epoch) / pd.Timedelta(days=1)
print(data.head())
print(f'tiene una longitud de {len(data)} filas')

# %%
ncfile = Dataset(f"{path}/{nombre_fichero}.nc", mode='w', format='NETCDF3_CLASSIC')

ncfile.title = nombre_fichero
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = n_dataset
ncfile.project = n_dataset; ncfile.source = 'In situ data collection'; ncfile.Conventions = 'CF-1.8'
ncfile.comment ='Coordinates for RAMBLA station(s) are approximate'

ncfile.createDimension('time', len(dias_desde_1970))
ncfile.createDimension('unit_char_len', max_length_param)
ncfile.createDimension('station_name', len(estaciones_np))
ncfile.createDimension('depth', len(depths_enteros))
ncfile.createDimension('depth_bottom', len(depths_bottom))

for dim in ncfile.dimensions.items():
    print(dim)

time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values  # Se asigna directamente

lat_var = ncfile.createVariable('latitude', np.float64,('station_name') )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var.grid_mapping = "crs"
lat_var[:] =  [37.8, 37.7833, 37.76667, 37.76667, 37.75, 37.75, 37.75, 
          37.73333 ,37.73333, 37.73333  , 37.71667 , 37.71667 , 37.71667 ,
            37.71667,37.7, 37.7,  37.7, 37.7, 37.68333, 37.68333, 
             37.68333, 37.66667,  37.66667, 37.65, 37.716111 ]

lon_var = ncfile.createVariable('longitude', np.float64,('station_name') )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] = [-0.78056, -0.78333, -0.8 , -0.76667, -0.81667 , -0.78333, -0.75, 
            -0.83333, -0.8, 0.76667, -0.85 , -0.81667, -0.78333,
           -0.75, -0.83333, -0.8, -0.78333, -0.75, -0.81667, -0.78333,
            -0.75, -0.8, -0.76667, -0.75, -0.86  ]

depth_var = ncfile.createVariable('depth',np.int8, ('depth',))
depth_var.units = 'meters'
depth_var.standard_name = 'depth'
depth_var.positive = 'down'
depth_var[:] = depths_enteros

depth_var = ncfile.createVariable('depth_bottom',np.int8, ('depth_bottom',))
depth_var.units = 'meters'
depth_var.long_name = 'depth at the bottom of the station'
depth_var.positive = 'down'
depth_var[:] = depths_bottom

value_var = ncfile.createVariable('nitrate', np.float32, ('time', 'station_name', 'depth'))
value_var.units = 'µmol/L'
value_var.standard_name = 'mole_concentration_of_nitrate_in_sea_water'
value_var.long_name = 'nitrate concentration in sea water'
value_var.missing_value = -9999
value_var.grid_mapping = "crs"
value_var.comment = ''

value_var2 = ncfile.createVariable('nitrate_at_bottom', np.float32, ('time', 'station_name'))
value_var2.units = 'µmol/L'
value_var2.standard_name = 'mole_concentration_of_nitrate_in_sea_water'
value_var2.long_name = 'Nitrate concentration in sea water'
value_var2.missing_value = -9999
value_var2.grid_mapping = "crs"
value_var2.comment = 'Nitrate at bottom of the station'

data_entera = data.loc[data['profundidad'].isin(depths_enteros)]


pivot = data.pivot_table(index='fecha',  columns=['estacion', 'profundidad'], values='valor')
# Todas las fechas
all_dates = pd.to_datetime(data.fecha.drop_duplicates()).sort_values()

# Reindexar pivot para que tenga todas las fechas
pivot = pivot.reindex(all_dates)

times = pd.to_datetime(data.fecha.drop_duplicates()).sort_values()
valor_array_3d = np.full((len(times), len(estaciones_np), len(depths_enteros)), np.nan)

for i, station in enumerate(estaciones):
    for j, depth in enumerate(depths_enteros):
        try:
            valor_array_3d[:, i, j] = pivot[(station, depth)].values
        except KeyError:
            pass  # si no hay datos

value_var[:,:,:] = valor_array_3d
# -----
data_fondo_mask = (
    data['profundidad'].isin(depths_bottom) |
    ((data['estacion'] == 'E14') & (data['profundidad'] == 5)) |
    ((data['estacion'] == 'E13') & (data['profundidad'] == 6))|
    ((data['estacion'] == 'E19') & (data['profundidad'] == 4))|
    ((data['estacion'] == 'RAMBLA') & (data['profundidad'] == 1))
)
# Aplicamos la máscara al DataFrame
data_fondo = data.loc[data_fondo_mask]
data_fondo = data_fondo.sort_values(by=["fecha", "profundidad", "estacion"]).reset_index(drop=True)
# Ahora sí podemos ordenar
valor_array_2d = np.full((len(times), len(estaciones)), np.nan)
for i, station in enumerate(estaciones):
    array_valores=[]
    for j, dia in enumerate(times):
        valor = data_fondo.loc[(data_fondo['estacion']==station) &
                                                 (data_fondo['fecha']==dia) ]['valor'].values
        if len(valor) > 0:
            array_valores.append(valor[0])
        else: 
            array_valores.append(np.nan)
    valor_array_2d[:, i] = array_valores

value_var2[:,:] = valor_array_2d
# -----

parameter_var = ncfile.createVariable('station', 'S1', ('station_name', 'unit_char_len'))
parameter_var.long_name = 'station'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)

bottom_depth_var = ncfile.createVariable('station_bottom_depth', np.float32, ('station_name',))
bottom_depth_var.units = 'meters'
bottom_depth_var.standard_name = 'sea_floor_depth_below_sea_surface'
bottom_depth_var.long_name = 'Depth of sea floor below sea surface at each station'
bottom_depth_var.positive = 'down'
bottom_depth_var[:] =[1.4, 5.1, 5.3, 5.3, 5.5, 6.2, 4.8, 5.1, 5.9, 6.1, 3.7, 5.3,6, 5, 3.5,  3.2, 5.7, 4.6, 4, 5.5, 5.5, 4.2, 5.3, 5.3, 1 ] # array de 24 valores, uno por estación


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


unit2 = dataset.variables['nitrate_at_bottom'][:]
unit = dataset.variables['nitrate'][:]
station = dataset.variables['station'][:]
tiempo = dataset.variables["time"][:]  # Días desde 1970

northing = dataset.variables['latitude'][:]
easting = dataset.variables['longitude'][:]
crs = dataset.variables['crs']
station_depth= dataset.variables['station_bottom_depth'][:]

print(station)
print(tiempo); print('-----------------')
print(unit); print('-----------------')
print(unit2); print('-----------------')

print(northing); print('-----------------')
print(easting); print('-----------------')
print(crs); print('-----------------')

fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
fechas = np.array(fechas) 
dataset.close()

# %%

generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')
# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Biogeochemical/nutrients/nitrate/'

shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')

# %%FIGURA
# Reemplazamos -9999 por NaN
unit_plot = unit.astype(float)

# Reemplazamos -9999 por NaN
unit_plot[unit_plot == -9999] = np.nan

for i, station in enumerate(estaciones):
    plt.figure(figsize=(12, 6))
    
    # Dibujar líneas para cada profundidad
    for j, depth in enumerate(depths_enteros):
        perfil = unit[:, i, j]
        if np.any(~np.isnan(perfil)):  # Solo dibuja si hay al menos un valor válido
            bootm_depth = station_depth[i]
            plt.plot(times, perfil, label=f'Depth {depth} m')

    
    # Dibujar línea de fondo
    bootm_depth = station_depth[i]
    plt.plot(times, unit2[:, i], 'k--', linewidth=2, label=f'Bottom {bootm_depth}')
    
    plt.title(f'Temporal Profile - Station {station}')
    plt.xlabel('Date')
    plt.ylabel('nitrate (µmol/l)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    file_name = f'{station}_temporal_profile.png'
    plt.savefig(os.path.join(path, file_name), dpi=300)
    plt.show()
    plt.close()
# %%
