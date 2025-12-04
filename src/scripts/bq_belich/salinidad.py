# %%
from netCDF4 import Dataset, stringtochar
import numpy as np
import pandas as pd
import shutil
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
from generar_txt import generar_txt
# %%
n_dataset ='BELICH'
nombre_fichero = 'BELICH_SAL'

data0 = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/scripts/bq_belich/Libro4_salinidad.xlsx',  dtype={  
"Valor": "float64", }, parse_dates= ['Fecha'],  usecols= ['Variable', 'Estación', 'Fecha','Valor','Profundidad' ])

print(data0.head())
print(f'tiene una longitud de {len(data0)} filas')
print(data0.tail())
data0 = data0.replace(np.nan, -9999) #reemplazo los nan por -9999 
# %%
max_length_param = len("A")
estaciones = ['A', 'B', 'C', 'M']
estaciones_np = np.array(estaciones, dtype=f'S{max(len(s) for s in estaciones)}')
depths = [0.5,2,4]

path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Hydrodynamics/salinity/BELICH/'
path_copia = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Hydrodynamics/salinity/'

data = data0.sort_values(by=["Fecha", "Estación"]).reset_index(drop=True)
print('la longitud antes es', len(data))
data = data.drop_duplicates()
print('la longitud despues es', len(data))
data['Fecha'] = pd.to_datetime(data['Fecha'], errors='coerce')
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (data.Fecha.drop_duplicates() - epoch) / pd.Timedelta(days=1)
print(data.head())
print(f'tiene una longitud de {len(data)} filas')
# %%
ncfile = Dataset(f"{path}/{nombre_fichero}.nc", mode='w', format='NETCDF3_CLASSIC')

ncfile.title = nombre_fichero
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = n_dataset
ncfile.project = 'BELICH'; ncfile.source = 'In situ data collection'; ncfile.Conventions = 'CF-1.8'

ncfile.createDimension('time', len(dias_desde_1970))
ncfile.createDimension('unit_char_len', max_length_param)
ncfile.createDimension('station_name', len(estaciones))
ncfile.createDimension('depth', len(depths))

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
lat_var[:] = [ 37.791433,  37.70949,  37.667317, 37.710278, ]

lon_var = ncfile.createVariable('longitude', np.float64,('station_name') )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] = [ -0.78155 , -0.7851 , -0.755215, -0.831111]

depth_var = ncfile.createVariable('depth',np.int8, ('depth',))
depth_var.units = 'meters'
depth_var.standard_name = 'depth'
depth_var.positive = 'down'
depth_var[:] = depths

value_var = ncfile.createVariable('seawater_salinity', np.float32, ('time', 'station_name', 'depth'))
value_var.units = '1'
value_var.standard_name = 'sea_water_practical_salinity'
value_var.long_name = 'Practical salinity of sea water'
value_var.missing_value = -9999
value_var.grid_mapping = "crs"
value_var.comment = 'salinity time series profiles collected at 4 stations identified by geographic coordinates (longitude, latitude). Reported in Practical Salinity Units (PSU), which are dimensionless.'

pivot = data.pivot_table(index='Fecha',  columns=['Estación', 'Profundidad'], values='Valor')
# Todas las fechas
all_dates = pd.to_datetime(data.Fecha.drop_duplicates()).sort_values()

# Reindexar pivot para que tenga todas las fechas
pivot = pivot.reindex(all_dates)

times = pd.to_datetime(data.Fecha.drop_duplicates()).sort_values()
valor_array_3d = np.full((len(times), len(estaciones), len(depths)), np.nan)

for i, station in enumerate(estaciones):
    for j, depth in enumerate(depths):
        try:
            valor_array_3d[:, i, j] = pivot[(station, depth)].values
        except KeyError:
            pass  # si no hay datos

value_var[:,:,:] = valor_array_3d

parameter_var = ncfile.createVariable('station', 'S1', ('station_name', 'unit_char_len'))
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


unit = dataset.variables['seawater_salinity'][:]
tiempo = dataset.variables["time"][:]  # Días desde 1970

northing = dataset.variables['latitude'][:]
easting = dataset.variables['longitude'][:]
crs = dataset.variables['crs']

print(tiempo); print('-----------------')
print(unit); print('-----------------')

print(northing); print('-----------------')
print(easting); print('-----------------')
print(crs); print('-----------------')

fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
fechas = np.array(fechas) 
dataset.close()

# %%

generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')
# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Hydrodynamics/salinity/'

shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')

# %%FIGURA
# Reemplazamos -9999 por NaN
unit_plot = unit.astype(float)

# Reemplazamos -9999 por NaN
unit_plot[unit_plot == -9999] = np.nan

# Definimos nombres y colores de estaciones
estaciones = ['A','B','C','M']
colors = {'A':'blue','B':'green','C':'orange','M':'red'}
profundidades = [0.5, 2, 4]

plt.figure(figsize=(12,6))

# Iteramos sobre estaciones y profundidades
for i_station, station in enumerate(estaciones):
    for i_depth, depth in enumerate(profundidades):
        plt.plot(fechas, unit_plot[:, i_station, i_depth], 
                 label=f'{station} {depth} m',
                 color=colors[station],
                 linestyle='-' if depth==0.5 else '--' if depth==2 else ':')

plt.xlabel('Fecha')
plt.ylabel('salinidad')
plt.title('Serie temporal de salinidad en el Mar Menor')
plt.legend(ncol=2, fontsize=9)
plt.grid(True)
plt.tight_layout()
plt.show()
# %%
