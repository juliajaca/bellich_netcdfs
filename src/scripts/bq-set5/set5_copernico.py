# %%
from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import pandas as pd
import shutil
import matplotlib.pyplot as plt
from geopy.distance import geodesic
import sys
sys.path.append("../")
from generar_txt import generar_txt
# %%
data = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/bq/set5/serie_chlacopernico.xlsx',  dtype={  
"Chla satelital corregida": "float64",}, parse_dates= ['FECHA'], ) #rows=3
data = data.sort_values(by='FECHA').reset_index(drop=True)
print(data.head())
print(f'tiene una longitud de {len(data)} filas')

# %%
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (data.FECHA - epoch) / pd.Timedelta(days=1)
# %%
nombre_fichero = 'IEO_SAT_CORR_CHL'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/chlorophyll/IEO_SAT_CORR/'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title = nombre_fichero
ncfile.institution = "Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'IEO_SAT_CORR' 
ncfile.project = 'Not associated with a specific project'
ncfile.source = 'Measurements from Copernicus satellite'
ncfile.Conventions = 'CF-1.8'

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

valores_con_nan = data['Chla satelital corregida']
valores_con_nan[np.isnan(valores_con_nan)] = -9999
print(valores_con_nan)

value_var = ncfile.createVariable('chlorophyll', np.float64, ('time',))
value_var.units= 'ug L-1'
value_var.standard_name = 'mass_concentration_of_chlorophyll_a_in_sea_water'
value_var.long_name = 'Chlorophyll-a Concentration in Sea Wate'
value_var.cell_methods= 'Mean chlorophyll concentration over a specified region from satellite. Methodology missing!'
value_var.missing_value = -9999
value_var.comment = 'BELA chlorophyll concentration algorithm used to estimate CHL. Corrected chlorophyll from satellite'
value_var[:] =valores_con_nan

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
value = dataset.variables["chlorophyll"][:] 
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
axes.set_ylabel('Clorofila ')
axes.set_title(f'Clorofila satelital corregida copernicus')
axes.grid(True)
axes.legend()

axes.set_xlabel('Fecha')
plt.tight_layout()
plt.show()

# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}.txt')
# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Biogeochemical/chlorophyll/IEO_SAT_CORR/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')
# %%