# %%
from netCDF4 import Dataset, stringtochar
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
from generar_txt import generar_txt
# %%
nombre_fichero = 'IEO_PRMAX'
data0 = pd.read_excel(
    'C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/atmosferico/set2/Precipitaciones_p99_max per year.xlsx',
    dtype={'E7020': 'float', 'E7025': 'float', 'E7026': 'float', 'E7029': 'float', 'Year': 'int'},
    usecols=['E7020', 'E7025', 'E7026', 'E7029', 'Year']
)

data1 = data0.sort_values(by=["Year"]).reset_index(drop=True)
data = data1.dropna(subset=data1.columns.difference(['Year']), how='all')
fechas = pd.to_datetime(data['Year'].astype(str) + '-01-01')

# Paso 2: Calcular días desde el epoch (1-1-1970)
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (fechas - epoch) / pd.Timedelta(days=1)

# %%
estaciones = data.columns[1:]
estaciones_np = np.array(estaciones, dtype=f'S{max(len(s) for s in estaciones)}')
max_length_param = max(len(s) for s in estaciones)

# %%
path= 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Atmospheric/precipitation/IEO_PR/'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')

ncfile.title=f'{nombre_fichero}'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'; ncfile.source = 'XXX'; ncfile.Conventions = 'CF-1.8'

ncfile.createDimension('time', len( data.Year.drop_duplicates() ))
ncfile.createDimension('unit_char_len', max_length_param)
ncfile.createDimension('station', len(estaciones))

for dim in ncfile.dimensions.items():
    print(dim)

time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
time_var[:] = dias_desde_1970.values  # Se asigna directamente

parameter_var = ncfile.createVariable('station_code', 'S1', ('station', 'unit_char_len'))
parameter_var.long_name = 'station'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)

value_var = ncfile.createVariable('precipitation', np.float64, ('time', 'station'))
value_var.standard_name = 'precipitation_flux'
value_var.units= 'XXX'
value_var[:, :] =  data[estaciones].to_numpy()   # o tambien data[estaciones].values

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

print(tiempo); print('-----------------')
print(unit); print('-----------------')
print(value); print('-----------------')
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
fechas = np.array(fechas) 
dataset.close()

# %%# %% PLOT
fig, axs = plt.subplots(2, 2, figsize=(16, 12), sharex=True, sharey=True)
axs = axs.flatten()

# Ajuste del rango Y (opcional: calcula el min y max global para todos los plots)
ymin, ymax = np.nanmin(value), np.nanmax(value)

for i, ax in enumerate(axs):
    ax.plot(fechas, value[:, i], label='precipitacion', color='blue', marker='s')
    ax.set_title(f"Estación {estaciones[i]}")
    ax.grid(True)
    ax.set_ylim(ymin, ymax)  # mismo rango en Y para todos

    # Mostrar ylabel solo en la primera columna
    if i % 2 == 0:
        ax.set_ylabel('precipitacion')
    else:
        ax.set_ylabel('')

    # Rotar las fechas del eje x solo si estamos en la última fila
    if i // 2 == 2:
        ax.tick_params(axis='x', rotation=90)

# Etiqueta general del eje X
fig.text(0.5, 0.04, 'Fecha', ha='center')

plt.tight_layout(rect=[0, 0.05, 1, 0.97])  # deja espacio para el eje X global
plt.savefig(f'{path}{nombre_fichero}.pdf', format='pdf')
plt.show()

# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')
# %%
