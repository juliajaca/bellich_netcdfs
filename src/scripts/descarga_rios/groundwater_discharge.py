# %%
from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# # %%
# lilia = Dataset('../datos/ejemplo_nc_lilia/set5_TEMP_IEO_SATELLITE.nc', "r")

data = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/rios/groundwater_discharge.xlsx',  )

# %%
path= 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/AportesContinentales/set1_runoff/'
ncfile = Dataset(f'{path}groundwater_discharge.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)
# Crear dimensiones
# Suponiendo que tienes 11 sectores con nombres de hasta 4 caracteres
sector_dim = ncfile.createDimension('sector', 11)
char_len_dim = ncfile.createDimension('char_len', 4)

print('---------')
print(ncfile)
for dim in ncfile.dimensions.items():
    print(dim)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='Groundwater Discharge'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon, but Np1 is discharging not in the lagoon but in Salinas de San Pedro, in the north'
ncfile.project = 'XXXX'
ncfile.source = 'Anexo 6. Resultado del modelo de flujo subterráneo para el año 2018/2019. Estudio para la cuantificaciónn de agua subterránea en el Mar Menor.'
ncfile.Conventions = 'CF-1.8'

# %%
# Define two variables with the same names as dimensions,
# a conventional way to define "coordinate variables".
sector_var = ncfile.createVariable('station_name', 'S1', ('sector','char_len'))
# sector_var.standard_name = 'region'
sector_var.long_name = 'Lagoon sector'

discharge_var = ncfile.createVariable('water_volume_transport', np.float64, ('sector',))
discharge_var.units = 'm3 d-1'
discharge_var.standard_name = "water_volume_transport_in_river_channel"
discharge_var.long_name = "Daily water discharge per lagoon sector"

# %% WRITING DATA 
max_length = data.Sector.str.len().max() #es 4
# sector_var[:, :] = np.array([list(s.ljust(4)) for s in data.Sector], dtype='S1') Esto es una forma

sector_var[:,:] = stringtochar(data.Sector.to_numpy(dtype=f'S{max_length}')) ## manual conversion to char array, https://unidata.github.io/netcdf4-python/
sector_var._Encoding = 'ascii' # this enables automatic conversion
print(sector_var[:])

# %%
discharge_var[:] = data.Discharge

# %%
ncfile.close(); print('Dataset is closed!')

# %%
dataset = Dataset(f'{path}groundwater_discharge.nc', "r")
print(dataset)
print(dataset.variables.keys())  # Ver las variables en el archivo
print('los discharge data son')
discharge_data = dataset.variables["water_volume_transport"][:]
print(discharge_data)
print('los sectores son ')
print(dataset.variables['station_name'][:])

print("\n🔹 Atributos de las Variables:")
for var_name in dataset.variables:
    print(f"\nVariable: {var_name}")
    for attr in dataset.variables[var_name].ncattrs():
        print(f"  {attr}: {dataset.variables[var_name].getncattr(attr)}")

print("\n🔹 Atributos Globales:")
for attr in dataset.ncattrs():
    print(f"{attr}: {dataset.getncattr(attr)}")

# %%
#  PLOT
# Leer los sectores y la variable discharge
sector_var = dataset.variables["station_name"][:]
sectores = [sector for sector in sector_var]

# sectores = ["".join(s.decode("utf-8") for s in sector) for sector in sector_var] #amtes del encoding
discharge_data = dataset.variables["water_volume_transport"][:]

# Cerrar el archivo
dataset.close()

plt.figure(figsize=(8, 5))
plt.bar(sectores, discharge_data, color="royalblue")

# Etiquetas
plt.xlabel("Sector")
plt.ylabel("Discharge (m³/d)")
plt.title("Flujo por sector")
plt.xticks(rotation=45)  # Rotar los nombres si están muy juntos
plt.grid(axis="y", linestyle="--", alpha=0.7)

# Mostrar el gráfico
plt.show()
# %%

# GENERAR TXT
ncfile = Dataset(f'{path}groundwater_discharge.nc', "r")
txt_filename = f"{path}groundwater_discharg_output.txt"

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
