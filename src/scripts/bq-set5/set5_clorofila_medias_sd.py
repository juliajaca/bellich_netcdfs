# %%
from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from geopy.distance import geodesic
import sys
sys.path.append("../")
from generar_txt import generar_txt
# %%
data0 = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/bq/set5/serieChla2_MEDIAS_SD.xlsx',  dtype={  
"Chla_in_situ": "float64",}, parse_dates= ['Date_in_situ'], ) #rows=3
print(data0.head())
print(f'tiene una longitud de {len(data0)} filas')

# %%
data = data0.dropna(subset=['Date_in_situ', 'Chla_in_situ'], how='all')
data = data.sort_values(by='Date_in_situ').reset_index(drop=True)
#  %%

epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (data.Date_in_situ - epoch) / pd.Timedelta(days=1)
# %%
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/set5/chl medias sd/'
ncfile = Dataset(f'{path}set5_clorofila_medias_sd.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='set 5 clorofila medias sd in situ'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'; ncfile.source = 'XXX'; ncfile.Conventions = 'CF-1.8'


# %%
# crear dimensiones
ncfile.createDimension('time', len(dias_desde_1970))

for dim in ncfile.dimensions.items():
    print(dim)

# %%
time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
time_var[:] = dias_desde_1970.values  # Se asigna directamente

value_var = ncfile.createVariable('chla_in_situ', np.float64, ('time',))
value_var.long_name = 'mass_concentration_of_chlorophyll_in_sea_water'
value_var.units= 'mg m-3'
value_var[:] = data['Chla_in_situ']

# %%
ncfile.close()

# %% COMPROBACION
dataset = Dataset(f'{path}set5_clorofila_medias_sd.nc', "r")
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
value = dataset.variables["chla_in_situ"][:] 
print(tiempo); print('-----------------')
print(value); print('-----------------')
# %%
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()
# %% PLOT
fechas = np.array(fechas)  # convertir a numpy array para indexar
orden = np.argsort(fechas)  # obtener índices ordenados

# Reordenar fechas y nitrato
fechas_ordenadas = fechas[orden]
v_ordenado = value[orden]

fig, axes = plt.subplots(nrows=1, figsize=(10, 10), sharex=True)
axes.plot(fechas_ordenadas, v_ordenado, marker='o')
axes.set_ylabel('Clorofila in situ ')
axes.set_title(f'Clorofila in situ medias sd')
axes.grid(True)
axes.legend()

axes.set_xlabel('Fecha')
plt.tight_layout()
plt.show()

# %%
generar_txt(f'{path}set5_clorofila_medias_sd.nc', f'{path}set5_clorofila_medias_sd_display.txt')
# %%
