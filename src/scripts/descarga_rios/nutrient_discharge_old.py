from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# %%
data = pd.read_csv('../datos/descarga_nutrientes/Albujon_nutrient_outlet_CHS.csv', sep=';', dtype={  
"Parameter": "object",
"Unit":"object",
"Value":'float64', }, parse_dates= ['Date'], )

for column in data.columns:
    print((data[column].iloc[0]))
    print(type(data[column].iloc[0]))
    print('-------------------------------')
# %%
#  conversión de fechas
# Fecha de referencia (1 de enero de 1970)
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (data.Date - epoch) / pd.Timedelta(days=1)
# Conversion de las µ
data["Unit"] = data["Unit"].str.replace("µ", "micro")
data['Unit'] = data['Unit'].str.replace('º','deg')
# %%
ncfile = Dataset('nutrient_discharge.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='Albujon nutrient discharge'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.poject = 'XXXX'
ncfile.source = 'Confederación Hidrográfica del Segura (CHS)'
ncfile.conventions = 'XXX'

# Crear dimensiones
max_length_param = data.Parameter.str.len().max() 
max_length_unit = data.Unit.str.len().max() 

# Obtener la longitud máxima de las unidades
ncfile.createDimension('unit_char_len', max_length_unit)
ncfile.createDimension('obs', len(data))  # Número total de mediciones (cada fila del CSV es una medición individual)
ncfile.createDimension('string_len', max_length_param)

for dim in ncfile.dimensions.items():
    print(dim)
# %%
time_var = ncfile.createVariable('time', np.float64, ('obs',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
time_var[:] = dias_desde_1970.values  # Se asigna directamente

parameter_var = ncfile.createVariable('parameter', 'S1', ('obs', 'string_len'))
parameter_var.standard_name = 'parameter'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(data.Parameter.to_numpy(dtype=f'S{max_length_param}'))

value_var = ncfile.createVariable('value', np.float64, ('obs',))
value_var.standard_name = 'Parameter value'
value_var[:] = data["Value"].values

unit_var = ncfile.createVariable('unit', 'S1', ('obs', 'unit_char_len'))
unit_var.standard_name = 'unit'
unit_var._Encoding = 'ascii'
unit_var[:,:] = stringtochar(data.Unit.to_numpy(dtype=f'S{max_length_unit}'))

# %%
ncfile.close()

# %% COMPROBACION
dataset = Dataset('nutrient_discharge.nc', "r")
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
parameter = dataset.variables["parameter"][:]    # 
unit = dataset.variables["unit"][:]    #  
value = dataset.variables["value"][:]    #


print(tiempo)
print('-----------------')
print(parameter)
print('-----------------')
print(unit)
print('-----------------')
print(value)
print('-----------------')

# %%
parametros_unicos = pd.Series(parameter).unique()

# %%
# Obtener lista única de parámetros
df = pd.DataFrame({"time": tiempo, "parameter": parameter, "unit": unit, "value": value})
# Verificar datos
print(df.head())
# %%
# Graficar cada parámetro en un subplot diferente
fig, axes = plt.subplots(len(parametros_unicos), 1, figsize=(10, 5 * len(parametros_unicos)), sharex=True)

for i, param in enumerate(parametros_unicos):
    ax = axes[i]
    sub_df = df.loc[df["parameter"] == param]
    ax.plot(pd.to_datetime(sub_df["time"], origin="1970-01-01", unit="D") , 
            sub_df["value"],
            marker="o",
            linestyle="-",
            label=f"{param} ({sub_df['unit'].iloc[0]})")
    
    ax.set_ylabel(f"{param} ({sub_df['unit'].iloc[0]})")
    ax.legend()
    ax.grid()

plt.xlabel("Fecha")
plt.tight_layout()
plt.show()
# %%
dataset.close()
# %%
ncfile = Dataset('nutrient_discharge.nc', "r")
txt_filename = "nutrient_discharge_output.txt"

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