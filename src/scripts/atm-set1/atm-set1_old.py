# %%
from netCDF4 import Dataset, stringtochar
import numpy as np
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

# %%

data = data.sort_values(by=["Fecha", "Hora"]).reset_index(drop=True)

data['FechaJunta'] = pd.to_datetime(data['Fecha'].astype(str) + ' ' + data['Hora'].astype(str), dayfirst=True,)
epoch = pd.Timestamp('1970-01-01')
n_dias_desde_1970 = data.FechaJunta.drop_duplicates() 

# %%
estaciones = data['Codest'].unique()
estaciones_np = np.array(estaciones, dtype=f'S{max(len(s) for s in estaciones)}')
max_length_param = max(len(s) for s in estaciones)

# %%
coordenadas = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/atmosferico/set1/Estac_coord.xlsx')
texto = "ED50 o a ETRS89? Gonzalo González Barberá 10/12/24 \n"
for i in range(len(coordenadas)):
    texto += f'{coordenadas.iloc[i].Estacion} {coordenadas.iloc[i].Nombre}: UTM X-->{coordenadas.iloc[i].UTM_X}, UTM Y-->{coordenadas.iloc[i].UTM_Y} \n'
# %%
ncfile = Dataset(f'{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')

ncfile.title=f'{nombre_fichero}'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'; ncfile.source = 'XXX'; ncfile.conventions = 'XXX'
ncfile.coordinates = texto


ncfile.createDimension('time', len(n_dias_desde_1970))
ncfile.createDimension('unit_char_len', max_length_param)
ncfile.createDimension('station', len(estaciones))

for dim in ncfile.dimensions.items():
    print(dim)

pivot = data.pivot_table(index='FechaJunta', columns='Codest', values='Prec')
valor_array = pivot.to_numpy()

value_var = ncfile.createVariable('precipitacion', np.float64, ('time', 'station'))
value_var.standard_name = 'precipitación'
value_var.unit= 'XXXX'
value_var[:,:] = valor_array

time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
dias_desde_1970 = pd.Series(pivot.index)
dias_desde_1970 = (dias_desde_1970 - epoch) / pd.Timedelta(days=1)
time_var[:] = dias_desde_1970.values  # Se asigna directamente

parameter_var = ncfile.createVariable('station', 'S1', ('station', 'unit_char_len'))
parameter_var.standard_name = 'station'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)


ncfile.close()

# %% COMPROBACION
dataset = Dataset(f'{nombre_fichero}.nc', "r")
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
unit = dataset.variables["station"][:]    #  
value = dataset.variables['precipitacion'][:] 

print(tiempo); print('-----------------')
print(unit); print('-----------------')
print(value); print('-----------------')
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
generar_txt(f'{nombre_fichero}.nc', f'{nombre_fichero}_display.txt')
# %%
