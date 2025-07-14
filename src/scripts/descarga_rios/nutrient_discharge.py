from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# %%
data = pd.read_csv('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/descarga_nutrientes/Albujon_nutrient_outlet_CHS.csv', sep=';', dtype={  
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
fechas_unicas = np.sort(data["Date"].unique())
fechas_unicas_ts = pd.to_datetime(fechas_unicas)

epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (fechas_unicas_ts - epoch) / pd.Timedelta(days=1)

# %%
path= 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/AportesContinentales/CHS_ALBUJON_V3/'
ncfile = Dataset(f'{path}CHS_ALBUJON_V3_COND.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='Albujon conductivity discharge'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'
ncfile.source = 'Confederación Hidrográfica del Segura (CHS)'
ncfile.Conventions = 'CF-1.8'

# Cuantos parameter hay
lista_params = (data["Parameter"].unique())
# 
ncfile.createDimension('time', len(dias_desde_1970))

for dim in ncfile.dimensions.items():
    print(dim)
# %%
pivot = data.pivot_table(index='Date', columns='Parameter', values='Value')
valor_array = pivot.to_numpy()

time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
time_var[:] = dias_desde_1970.values  # Se asigna directamente

conductivity_var = ncfile.createVariable('conductivity', np.float64, ('time',))
conductivity_var.long_name = 'sea water electrical conductivity at 20 ºC'
conductivity_var.units= 'uS cm-1'
conductivity_var[:] = valor_array[:,0]

# %%
ncfile.close()

# %% COMPROBACION
dataset = Dataset(f'{path}CHS_ALBUJON_V3_COND.nc', "r")
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
# nitrate = dataset.variables["nitrate"][:]    # 
# fofate = dataset.variables["fosfate"][:]    #  
# conductivity = dataset.variables["conductivity"][:]    #
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
fechas = np.array(fechas) 

# print(tiempo[0:10])
# print('-----------------')
# print(nitrate[0:10])


# # %%
# Graficar cada parámetro en un subplot diferente
fig, axes = plt.subplots(3, 1, figsize=(10, 5 * (3)), sharex=True)
nombres = ['nitrato', 'fosfato', 'conductividad']
for i, param in enumerate([valor_array[:,0], valor_array[:,1], valor_array[:,2]]):
    ax = axes[i]
    ax.plot(fechas , 
            param,
            marker="o",
            linestyle="-")
    
    ax.set_ylabel(f"{nombres[i]} ")
    ax.legend()
    ax.grid()

plt.xlabel("Fecha")
plt.tight_layout()
plt.show()
# %%
dataset.close()
# %%
ncfile = Dataset(f'{path}CHS_ALBUJON_V3_COND.nc', "r")
txt_filename = f"{path}CHS_ALBUJON_V3_COND.txt"

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
"""
>==>    >=> >=> >===>>=====> >======>           >>       >===>>=====> >=======>
>> >=>  >=> >=>      >=>     >=>    >=>        >>=>           >=>     >=>
>=> >=> >=> >=>      >=>     >=>    >=>       >> >=>          >=>     >=>
>=>  >=>>=> >=>      >=>     >> >==>         >=>  >=>         >=>     >=====>
>=>   > >=> >=>      >=>     >=>  >=>       >=====>>=>        >=>     >=>
>=>    >>=> >=>      >=>     >=>    >=>    >=>      >=>       >=>     >=>
>=>     >=> >=>      >=>     >=>      >=> >=>        >=>      >=>     >=======>

"""
path= 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/AportesContinentales/CHS_ALBUJON_V3/'
ncfile = Dataset(f'{path}CHS_ALBUJON_V3_NO2.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='Albujon nitrate discharge'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'
ncfile.source = 'Confederación Hidrográfica del Segura (CHS)'
ncfile.Conventions = 'CF-1.8'

# Cuantos parameter hay
lista_params = (data["Parameter"].unique())
# 
ncfile.createDimension('time', len(dias_desde_1970))

for dim in ncfile.dimensions.items():
    print(dim)
# %%
pivot = data.pivot_table(index='Date', columns='Parameter', values='Value')
valor_array = pivot.to_numpy()

time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
time_var[:] = dias_desde_1970.values  # Se asigna directamente

nitrate_var = ncfile.createVariable('nitrate', np.float64, ('time',))
nitrate_var.long_name = 'mass concentration of nitrate in sea water'
nitrate_var.units= 'mg L-1'
nitrate_var[:] = valor_array[:,1]

# %%
ncfile.close()

# %% COMPROBACION
dataset = Dataset(f'{path}CHS_ALBUJON_V3_NO2.nc', "r")
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
dataset.close()
# %%
ncfile = Dataset(f'{path}CHS_ALBUJON_V3_NO2.nc', "r")
txt_filename = f"{path}CHS_ALBUJON_V3_NO2.txt"

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

"""
>=======>     >===>        >=>>=>   >=======>       >>       >===>>=====>     >===>
>=>         >=>    >=>   >=>    >=> >=>            >>=>           >=>       >=>    >=>
>=>       >=>        >=>  >=>       >=>           >> >=>          >=>     >=>        >=>
>=====>   >=>        >=>    >=>     >=====>      >=>  >=>         >=>     >=>        >=>
>=>       >=>        >=>       >=>  >=>         >=====>>=>        >=>     >=>        >=>
>=>         >=>     >=>  >=>    >=> >=>        >=>      >=>       >=>       >=>     >=>
>=>           >===>        >=>>=>   >=>       >=>        >=>      >=>         >===>

"""
path= 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/AportesContinentales/CHS_ALBUJON_V3/'
ncfile = Dataset(f'{path}CHS_ALBUJON_V3_PO4.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='Albujon phosphate discharge'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'
ncfile.source = 'Confederación Hidrográfica del Segura (CHS)'
ncfile.Conventions = 'CF-1.8'

# Cuantos parameter hay
lista_params = (data["Parameter"].unique())
# 
ncfile.createDimension('time', len(dias_desde_1970))

for dim in ncfile.dimensions.items():
    print(dim)
# %%
pivot = data.pivot_table(index='Date', columns='Parameter', values='Value')
valor_array = pivot.to_numpy()

time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
time_var[:] = dias_desde_1970.values  # Se asigna directamente


fosfate_var = ncfile.createVariable('phosphate', np.float64, ('time',))
fosfate_var.long_name = 'mass concentration of phosphate in sea water'
fosfate_var.units= 'mg L-1'
fosfate_var[:] = valor_array[:,2]
# %%
ncfile.close()

# %% COMPROBACION
dataset = Dataset(f'{path}CHS_ALBUJON_V3_PO4.nc', "r")
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


# %%
dataset.close()
# %%
ncfile = Dataset(f'{path}CHS_ALBUJON_V3_PO4.nc', "r")
txt_filename = f"{path}CHS_ALBUJON_V3_PO4.txt"

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
