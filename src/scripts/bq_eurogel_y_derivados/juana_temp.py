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
n_dataset ='IEO_COMU_2'
nombre_fichero = 'IEO_COMU_2_TEMP_SW' #TESIS DE JUANA CANO

data0 = pd.read_excel('C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Physico-chemical/IEO_EUROGEL/060325_resumen_historicos.xlsx',  dtype={  
"valor": "float64", }, parse_dates= ['fecha'],  usecols= ['var', 'estacion', 'fecha','profundidad','valor' , 'muestreo', 'latitud', 'longitud'])

# %%
oxy = data0.loc[(data0['var'] == 'temp') & (data0['muestreo'] == 'Juana Cano')]
print(data0.head())
print(f'tiene una longitud de {len(data0)} filas')
oxy = oxy.replace(np.nan, -9999) #reemplazo los nan por -9999 
print(oxy.tail())
# %%
max_length_param = len("E1")
estaciones = ['E1',  'E2', 'E3', 'E4', ]
estaciones_np = np.array(estaciones, dtype=f'S{max(len(s) for s in estaciones)}')


# %%
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Hydrodynamics/temperature/IEO_COMU_2/'
path_copia = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Hydrodynamics/temperature/'

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
ncfile.institution= "Oceanographic Center of Murcia (COMU), Spain"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = n_dataset
ncfile.project ='OSTRAS'
ncfile.source = 'In situ data collection'
ncfile.Conventions = 'CF-1.8'

ncfile.createDimension('time', len(dias_desde_1970))
ncfile.createDimension('unit_char_len', max_length_param)
ncfile.createDimension('station_name', len(estaciones_np))
ncfile.createDimension('depth', len([0,5.25]))


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
lat_var[:] =  [ 37.79604, 37.75158, 37.67723, 37.70578,]

lon_var = ncfile.createVariable('longitude', np.float64,('station_name') )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] = [ -0.781283,-0.784307, -0.786973,-0.824049]

depth_var = ncfile.createVariable('depth',np.float32, ('depth',))
depth_var.units = 'meters'
depth_var.standard_name = 'depth'
depth_var.positive = 'down'
depth_var.comment = 'Depth values indicate position only (not measured depths): 0 m for surface, 5.25 m for seabed.'
depth_var[:] = [0, 5.25]

value_var = ncfile.createVariable('seawater_temperature', np.float32, ('time', 'station_name','depth'))
value_var.units = 'degree_Celsius'
value_var.standard_name = 'sea_water_temperature'
value_var.long_name= 'Sea water temperature'
value_var.missing_value = -9999
value_var.grid_mapping = "crs"
value_var.cell_methods  = 'time: mean'
value_var.comment = ''

pivot = data.pivot_table(index='fecha',  columns=['estacion', 'profundidad'], values='valor')
# Todas las fechas
all_dates = pd.to_datetime(data.fecha.drop_duplicates()).sort_values()

# Reindexar pivot para que tenga todas las fechas
pivot = pivot.reindex(all_dates)

times = pd.to_datetime(data.fecha.drop_duplicates()).sort_values()
valor_array_3d = np.full((len(times), len(estaciones), len([0,1])), np.nan)

for i, station in enumerate(estaciones):
    print(station)
    for j, depth in enumerate(['superficie','fondo']):
       
        try:
            valor_array_3d[:, i, j] = pivot[(station, depth)].values
        except KeyError:
            pass  # si no hay datos

value_var[:,:,:] = valor_array_3d
# -----
# %%
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

unit = dataset.variables['seawater_temperature'][:]
stations = dataset.variables['station'][:]
tiempo = dataset.variables["time"][:]  # Días desde 1970
prof = dataset.variables['depth'][:]

northing = dataset.variables['latitude'][:]
easting = dataset.variables['longitude'][:]
crs = dataset.variables['crs']

print(stations)
print(tiempo); print('-----------------')
print(unit); print('-----------------')

print(northing); print('-----------------')
print(easting); print('-----------------')
print(prof); print('-----------------')
print(crs); print('-----------------')

fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
fechas = np.array(fechas) 
dataset.close()

# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')
# %%
shutil.copy(f'{path}{nombre_fichero}.nc',f'{path_copia}{nombre_fichero}.nc')

# %%FIGURA
# Reemplazamos -9999 por NaN
unit_plot = unit.astype(float)

# Reemplazamos -9999 por NaN
unit_plot[unit_plot == -9999] = np.nan

for i_station, station in enumerate(stations):

    plt.figure(figsize=(12,6))   # ← figura nueva para cada estación

    # Para cada profundidad, trazamos la serie de esa estación
    for i_depth, depth in enumerate(prof):
        plt.plot(fechas,
                 unit_plot[:, i_station, i_depth],
                 label=f'{depth} m')

    plt.xlabel('Fecha')
    plt.ylabel('temperatura ')
    plt.title(f'Serie temporal de temperatura en {station} en {northing[i_station]} y {easting[i_station]}')
    plt.legend(title="Profundidad", fontsize=9)
    plt.grid(True)
    plt.tight_layout()

    # Nombre único por estación
    file_name = f'{station}_temporal_profile.png'
    plt.savefig(os.path.join(path, file_name), dpi=300)

    plt.show()
    plt.close()
# %%
