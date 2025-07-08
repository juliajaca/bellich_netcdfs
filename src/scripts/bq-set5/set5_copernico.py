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
data = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/bq/set5/serie_chlacopernico.xlsx',  dtype={  
"Chla satelital corregida": "float64",}, parse_dates= ['FECHA'], ) #rows=3
data = data.sort_values(by='FECHA').reset_index(drop=True)
print(data.head())
print(f'tiene una longitud de {len(data)} filas')

# %%
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (data.FECHA - epoch) / pd.Timedelta(days=1)
# %%
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/chlorophyll/IEO_SAT_CORR/'
ncfile = Dataset(f'{path}set5_clorofila_copernicus2.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='Chlorophyll BELA measurements'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'IEO_SAT_CORR'; ncfile.source = 'Measurements from Copernicus satellite'; ncfile.Conventions = 'CF-1.8'


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

value_var = ncfile.createVariable('chlorophyll', np.float64, ('time',))
value_var.standard_name = 'chlorophyll_concentration_in_sea_water'
value_var.long_name = 'Area-averaged chlorophyll concentration from satellite'
value_var.units= 'mg m-3'
value_var.cell_methods= 'Mean chlorophyll concentration over a specified region from satellite. Methodology missing!'
value_var.comment = 'BELA chlorophyll concentration algorithm used to estimate CHL. Corrected chlorophyll from satellite'
value_var[:] = data['Chla satelital corregida']

# %%
ncfile.close()

# %% COMPROBACION
dataset = Dataset(f'{path}set5_clorofila_copernicus2.nc', "r")
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
generar_txt(f'{path}set5_clorofila_copernicus2.nc', f'{path}set5_clorofila_copernicus.txt')
# %%
