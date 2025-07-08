# %%
from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# %%
data = pd.read_csv('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/rios/Albujon_flow_hourly_CHS.csv', sep =',', dtype={  
"Flow_m3_s": "float64",}, parse_dates= ['Date_time'])
# %%
print(type(data.Date_time.iloc[0]))
print(type(data.Flow_m3_s.iloc[1]))
# %%
# Fecha de referencia (1 de enero de 1970)
epoch = pd.Timestamp('1970-01-01')

# Calcular la diferencia en días
dias_desde_1970 = (data.Date_time - epoch) / pd.Timedelta(days=1)

# %%
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/AportesContinentales/set3_runoff/'
ncfile = Dataset(f'{path}albujon_flow_hourly.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)
# Crear dimensiones
datetime_dim = ncfile.createDimension('time', len(data))
# flow_dim = ncfile.createDimension('flow', len(data))
print('---------')
print(ncfile)
for dim in ncfile.dimensions.items():
    print(dim)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='Albujon flow hourly CHS'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'
ncfile.source = 'Confederación Hidrográfica del Segura (CHS)'
ncfile.Conventions = 'CF-1.8'

# %%
datetime_var = ncfile.createVariable('time', np.float64, ('time'))
datetime_var.units= "days since 1970-01-01 00:00:0"
datetime_var.standard_name = "time"
datetime_var.calendar = 'gregorian'
datetime_var[:] = dias_desde_1970.values

flow_var = ncfile.createVariable('water_volume_transport', np.float64, ('time',))
flow_var.units = 'm3 s-1'
flow_var.standard_name = 'water_volume_transport_in_river_channel'
flow_var[:] = data["Flow_m3_s"].values
ncfile.close()

# %%
#  COMPROBACION
dataset = Dataset(f'{path}albujon_flow_hourly.nc', "r")
print(dataset)
print(dataset.variables.keys())  # Ver las variables en el archivo
print('los fechas data son')
discharge_data = dataset.variables["time"][:]
print(discharge_data)
print('los flows son ')
print(dataset.variables['water_volume_transport'][:])

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
flujo = dataset.variables["water_volume_transport"][:]    # Flujo en m3/s

# Convertir tiempo a formato datetime para mejor visualización
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")

# Cerrar el archivo
dataset.close()

# Hacer el plot
plt.figure(figsize=(10, 5))
plt.plot(fechas, flujo, marker="o", linestyle="-", color="b", label="Flujo (m³/s)")
plt.xlabel("Fecha")
plt.ylabel("Flujo (m³/s)")
plt.title("Serie Temporal del Flujo")
plt.xticks(rotation=45)
plt.legend()
plt.grid()
plt.show()
# %%
# GENERAR TXT
ncfile = Dataset(f'{path}albujon_flow_hourly.nc', "r")
txt_filename = f"{path}albujon_flow_hourly_output.txt"

with open(txt_filename, "w") as f:
    # Escribir el formato
    f.write(f"Format:\n\t{ncfile.file_format}\n")

    # Escribir atributos globales
    f.write("\nGlobal Attributes:\n")
    for attr in ncfile.ncattrs():
        f.write(f"\t{attr:<15} = '{getattr(ncfile, attr)}'\n")

    # Escribir dimensiones
    f.write("\nDimensions:\n")
    for dim_name, dim in ncfile.dimensions.items():
        f.write(f"\t{dim_name:<15} = {len(dim)}\n")

    # Escribir variables
    f.write("\nVariables:\n")
    for var_name, var in ncfile.variables.items():
        # Obtener dimensiones y tipo de datos
        dimensions = ", ".join(var.dimensions)
        dtype = var.dtype

        # Escribir información de la variable
        f.write(f"    {var_name:<15}\n")
        f.write(f"\tSize: {var.shape}\n")
        f.write(f"\tDimensions: {dimensions}\n")
        f.write(f"\tDatatype: {dtype}\n")
        
        # Escribir atributos de la variable
        if var.ncattrs():
            f.write("\tAttributes:\n")
            for attr in var.ncattrs():
                f.write(f"\t\t{attr:<15} = '{getattr(var, attr)}'\n")

# Cerrar el archivo NetCDF
ncfile.close()

print(f"Archivo '{txt_filename}' generado con éxito.")
# %%