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
nombre_fichero = 'precipitacion_SIAM_GGB'
estA = pd.read_csv('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/atmosferico/set1/Precip_SIAM_16_24_MMa.csv',  decimal=',', sep=';', parse_dates= ['Fecha'],) #son todo strings

estB = pd.read_csv('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/atmosferico/set1/Precip_SIAM_16_24_MMb.csv',  decimal=',', sep=';', parse_dates= ['Fecha'],)
data = pd.concat([estA, estB], ignore_index=True)

# orden_deseado = data["Codest"].unique()
# %%

data = data.sort_values(by=["Fecha", "Hora"]).reset_index(drop=True)

data['FechaJunta'] = pd.to_datetime(data['Fecha'].astype(str) + ' ' + data['Hora'].astype(str), dayfirst=True,)
epoch = pd.Timestamp('1970-01-01')
n_dias_desde_1970 = data.FechaJunta.drop_duplicates() 


# %%
estaciones = data['Codest'].unique()
estaciones_np = np.array(estaciones, dtype=f'S{max(len(s) for s in estaciones)}')
max_length_param = max(len(s) for s in estaciones)

coordenadas = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/atmosferico/set1/Estac_coord.xlsx')
coordenadas['Estacion'] = pd.Categorical(coordenadas['Estacion'], categories=estaciones, ordered=True)
coordenadas = coordenadas.sort_values('Estacion').reset_index(drop=True)

# %%
nombre_fichero = 'SIAM_PR_MEAN_HOURLY'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Atmospheric/precipitation/SIAM/'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')

ncfile.title= nombre_fichero
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'SIAM'
ncfile.project = 'Not associated with a specific project'
ncfile.source = 'In situ data collection'
ncfile.Conventions = 'CF-1.8'

'''
SIAM Servicio agrometeorológico de Murcia en la cuenca vertiente del Mar Menor, 12 stations, hourly sampling. (dataset: Precip_SIAM_16_24_Mma.xls, Precip_SIAM_16_24_Mmb.xls,  Estac_coord.xls, responsible: Gonzalo González Barberá, period covered: 1/01/2016 to 30/11/2024 pero varia con la estacion, path: Datos_MM_Art_2025\Atmosfericas\GGB, variables measured: hourly precipitation).

dataset_id = SIAM
project = Not associated with a specific project
'''

ncfile.createDimension('time', len(n_dias_desde_1970))
ncfile.createDimension('unit_char_len', max_length_param)
ncfile.createDimension('station', len(estaciones))
nombres_completos = [coordenadas.iloc[i].Nombre for i in range(len(coordenadas))]
max_len_nombre = max(len(s) for s in nombres_completos)
ncfile.createDimension('name_full_strlen', max_len_nombre)

for dim in ncfile.dimensions.items():
    print(dim)

pivot = data.pivot_table(index='FechaJunta', columns='Codest', values='Prec')
valor_array = pivot.to_numpy()

time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
dias_desde_1970 = pd.Series(pivot.index)
dias_desde_1970 = (dias_desde_1970 - epoch) / pd.Timedelta(days=1)
time_var[:] = dias_desde_1970.values  # Se asigna directamente

# Latitud y longitud
lat_var = ncfile.createVariable('station_northing', np.float64, ('station',))
lat_var.units = 'm'
lat_var.long_name = 'UTM northing'
lat_var.grid_mapping = "crs"
lat_var[:] =  [coordenadas.iloc[i].UTM_Y for i in range(len(coordenadas))]

lon_var = ncfile.createVariable('station_easting', np.float64, ('station',))
lon_var.units = 'm'
lon_var.long_name = 'UTM easting'
lon_var.grid_mapping = "crs"
lon_var[:] =  [coordenadas.iloc[i].UTM_X for i in range(len(coordenadas))]

valores_con_nan = valor_array
valores_con_nan[np.isnan(valores_con_nan)] = -9999
value_var = ncfile.createVariable('precipitation', np.float64, ('time', 'station'))
value_var.units= 'mm'
value_var.standard_name = 'lwe_thickness_of_precipitation_amount'
value_var.long_name = 'Annual maximum precipitation (liquid water equivalent)'
value_var.cell_methods= 'time: mean'
value_var.missing_value = -9999
value_var.grid_mapping = "crs"
value_var[:,:] = valores_con_nan

parameter_var = ncfile.createVariable('station_code', 'S1', ('station', 'unit_char_len'))
parameter_var.long_name = 'station code'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)

nombre_var = ncfile.createVariable('station_name', 'S1', ('station', 'name_full_strlen'))
nombre_var.long_name = 'station name and location'
nombre_var._Encoding = 'ascii'
nombre_var[:] = stringtochar(np.array(nombres_completos, dtype=f'S{max_len_nombre}'))

crs = ncfile.createVariable('crs', 'i')
crs.grid_mapping_name = "transverse_mercator"
crs.projection = "UTM"
crs.long_name = "ETRS89 / UTM zone 30N"
crs.epsg_code = "EPSG:25830"
crs.semi_major_axis = 6378137.0
crs.inverse_flattening = 298.257222101
crs.comment = "The UTM coordinates might refer to either ED50 or ETRS89; ETRS89 (EPSG:25830) has been assumed by default"

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
unit = dataset.variables["station_code"][:]    #  
value = dataset.variables['precipitation'][:] 
station = dataset.variables['station_name'][:]

print(tiempo); print('-----------------')
print(unit); print('-----------------')
print(value[0:10]); print('-----------------')
print(station); print('-----------------')
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
fechas = np.array(fechas) 
dataset.close()

# %% PLOT
fig, axs = plt.subplots(4, 3, figsize=(16, 12), sharex=True, sharey=True)
axs = axs.flatten()

# Ajuste del rango Y (opcional: calcula el min y max global para todos los plots)
ymin, ymax = np.nanmin(value), np.nanmax(value)

for i, ax in enumerate(axs):
    ax.plot(fechas, value[:, i], label='precipitacion', color='blue', marker='s')
    ax.set_title(f"Estación {estaciones[i]}")
    ax.grid(True)
    ax.set_ylim(ymin, ymax)  # mismo rango en Y para todos

    # Mostrar ylabel solo en la primera columna
    if i % 3 == 0:
        ax.set_ylabel('precipitacion')
    else:
        ax.set_ylabel('')

    # Rotar las fechas del eje x solo si estamos en la última fila
    if i // 3 == 3:
        ax.tick_params(axis='x', rotation=90)

# Etiqueta general del eje X
fig.text(0.5, 0.04, 'Fecha', ha='center')

plt.tight_layout(rect=[0, 0.05, 1, 0.97])  # deja espacio para el eje X global
plt.savefig(f'{nombre_fichero}.pdf', format='pdf')
plt.show()

# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')
# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Atmospheric/precipitation/SIAM/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')
# %%