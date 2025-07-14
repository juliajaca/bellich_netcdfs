# %%
from netCDF4 import Dataset, stringtochar
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
from generar_txt import generar_txt
# %% cargar columnas
cols = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/atmosferico/set3/Radiacion 1996_2022.xlsx', nrows=0).columns
columnas_que_quiero =cols[list(range(1, 4)) + list(range(17, len(cols)))]
dtypes_dict = {col: "int" if col in ['AÑO','MES', 'DIA'] else 'float64' 
               for col in columnas_que_quiero}
# %%
# nombre_fichero = 'radiacion_1996_2002'
nombre_fichero = 'AEMET_SRAD'

data0 = pd.read_excel(
    'C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/atmosferico/set3/Radiacion 1996_2022.xlsx',
    dtype= dtypes_dict,
    usecols= columnas_que_quiero
)
# %%
estaciones = data0.columns[3:]
estaciones_np = np.array(estaciones, dtype=f'S{max(len(s) for s in estaciones)}')
max_length_param = max(len(s) for s in estaciones)

data0['Fecha'] = pd.to_datetime(data0['AÑO'].astype(str) + '-' + data0['MES'].astype(str) + '-' + data0['DIA'].astype(str),
    format='%Y-%m-%d')

data = data0.sort_values(by=["Fecha"]).reset_index(drop=True)
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (data.Fecha - epoch) / pd.Timedelta(days=1)

# %%
# path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Atmospheric/solarRadiation/set3-julia/'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Atmospheric/solarRadiation/AEMET/'

ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')

ncfile.title=f'{nombre_fichero}'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'; ncfile.source = 'XXX'; ncfile.Conventions = 'CF-1.8'
ncfile.indicativo = '7178I'
ncfile.altitud = '61'
ncfile.ind_syn = '8430'
ncfile.nom_eg = 'AEMET'

ncfile.createDimension('time', len( dias_desde_1970 ))
ncfile.createDimension('unit_char_len', max_length_param)
# ncfile.createDimension('variable_radiacion', len(estaciones))

for dim in ncfile.dimensions.items():
    print(dim)

time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
time_var[:] = dias_desde_1970.values  # Se asigna directamente

for column in data[estaciones]:
    print(column)
    value_var = ncfile.createVariable(column, np.float64, ('time',))
    value_var.long_name = column
    value_var.units= 'XXXX'
    value_var[:] = data[estaciones][column]

lat_var = ncfile.createVariable('northing', np.float64, )
lat_var.units = 'm'
lat_var.long_name = 'UTM northing'
lat_var.grid_mapping = "crs"
lat_var[:] =  4207610

lon_var = ncfile.createVariable('easting', np.float64, )
lon_var.units = 'm'
lon_var.long_name = 'UTM easting'
lon_var.grid_mapping = "crs"
lon_var[:] =  660598

crs = ncfile.createVariable('crs', 'i')
crs.grid_mapping_name = "transverse_mercator"
crs.projection = "UTM"
crs.long_name = "ETRS89 / UTM zone 30N"
crs.epsg_code = "EPSG:25830"

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
unit = dataset.variables["RDIR05"][:]    #  
unit2 = dataset.variables["RDIF06"][:]    #  
unit3 = dataset.variables['PTJERGL'][:] 
northing = dataset.variables['northing'][:]
easting = dataset.variables['easting'][:]
crs = dataset.variables['crs']

print(tiempo); print('-----------------')
print(unit); print('-----------------')
print(unit2); print('-----------------')
print(unit3); print('-----------------')
print(northing); print('-----------------')
print(easting); print('-----------------')
print(crs); print('-----------------')

fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
fechas = np.array(fechas) 
dataset.close()

# %% PLOT
# estaciones = np.array([
#     'RDIR05', 'RDIR06', 'RDIR07', 'RDIR08', 'RDIR09', 'RDIR10', 'RDIR11', 'RDIR12',
#     'RDIR13', 'RDIR14', 'RDIR15', 'RDIR16', 'RDIR17', 'RDIR18', 'RDIR19', 'RDIR20',
#     'RDIF05', 'RDIF06', 'RDIF07', 'RDIF08', 'RDIF09', 'RDIF10', 'RDIF11', 'RDIF12',
#     'RDIF13', 'RDIF14', 'RDIF15', 'RDIF16', 'RDIF17', 'RDIF18', 'RDIF19', 'RDIF20',
#     'RGLO05', 'RGLO06', 'RGLO07', 'RGLO08', 'RGLO09', 'RGLO10', 'RGLO11', 'RGLO12',
#     'RGLO13', 'RGLO14', 'RGLO15', 'RGLO16', 'RGLO17', 'RGLO18', 'RGLO19', 'RGLO20',
#     'RDIRDIA', 'RDIFDIA', 'RGLODIA', 'PTJERGL', 'PTJERGF'
# ])


# plots_por_pagina = 4
# filas, columnas = 2, 2
# total_series = len(estaciones)
# total_paginas = math.ceil(total_series / plots_por_pagina)

# for pagina in range(total_paginas):
#     fig, axs = plt.subplots(filas, columnas, figsize=(18, 10), sharex=True)
#     axs = axs.flatten()

#     start_idx = pagina * plots_por_pagina
#     end_idx = min(start_idx + plots_por_pagina, total_series)

#     for i, idx in enumerate(range(start_idx, end_idx)):
#         ax = axs[i]
#         ax.plot(fechas, value[:, idx], label=estaciones[idx], marker='o', linewidth=1)
#         ax.set_title(estaciones[idx])
#         ax.grid(True)

#         if i % columnas == 0:
#             ax.set_ylabel(f'Radiacion')
#         else:
#             ax.set_ylabel('')

#         if i // columnas == filas - 1:
#             ax.tick_params(axis='x', rotation=90)

#         ax.legend(fontsize='small', loc='upper right')

#     # Ocultar ejes vacíos si quedan huecos
#     for j in range(end_idx - start_idx, len(axs)):
#         fig.delaxes(axs[j])

#     fig.text(0.5, 0.04, 'Fecha', ha='center')
#     fig.suptitle(f'Página {pagina + 1}', fontsize=16)
#     plt.tight_layout(rect=[0, 0.05, 1, 0.95])
#     plt.savefig(f'{path}{nombre_fichero}_pagina_{pagina + 1}.pdf', format='pdf')
#     plt.close()
# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')
# %%
