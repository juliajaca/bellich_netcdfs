from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import pandas as pd
import shutil
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
nombre_fichero='CHS_ALBUJON_V3_COND'
path= 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Runoff/CHS_ALBUJON_V3/'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title=f'{nombre_fichero}'
ncfile.institution="Confederación Hidrográfica del Segura (CHS)"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'CHS_ALBUJON_V3'
ncfile.project = 'Not associated with a specific project'
ncfile.source = 'In situ data collection'
ncfile.Conventions = "CF-1.8"

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
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values  # Se asigna directamente

# Lo paso a mS cm --> 1 mS/cm= 1000 μS/cm
valores_con_nan = valor_array[:,0] / 1000
valores_con_nan[np.isnan(valores_con_nan)] = -9999 
# %%
conductivity_var = ncfile.createVariable('conductivity', np.float64, ('time',))
# conductivity_var.units= 'uS cm-1' # unidades anteriores, no son estandar
conductivity_var.units = 'mS cm-1'
conductivity_var.standard_name = "water_body_electrical_conductivity"
conductivity_var.long_name = 'Water body electrical conductivity'
conductivity_var.cell_methods= "time: mean" # no se si es mean o puntual
conductivity_var.missing_value = -9999 
conductivity_var.comment = 'Measured in surface runoff from Rambla del Albujón, a freshwater input to the Mar Menor'
conductivity_var[:] = valores_con_nan

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
ncfile = Dataset(f'{path}{nombre_fichero}.nc', "r")
txt_filename = f"{path}{nombre_fichero}.txt"

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
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Runoff/CHS_ALBUJON_V3/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')
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
nombre_fichero = 'CHS_ALBUJON_V3_NO2'
path= 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Runoff/CHS_ALBUJON_V3/'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title=f'{nombre_fichero}'
ncfile.institution="Confederación Hidrográfica del Segura (CHS)"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'CHS_ALBUJON_V3'
ncfile.project = 'Not associated with a specific project'
ncfile.source = 'In situ data collection'
ncfile.Conventions = "CF-1.8"

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
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values  # Se asigna directamente

valores_con_nan = (valor_array[:,1] *1000) / 62.0049 #paso de mg L-1 a umol L-1 
valores_con_nan[np.isnan(valores_con_nan)] = -9999

nitrate_var = ncfile.createVariable('nitrate', np.float64, ('time',))
nitrate_var.units= 'umol L-1'
nitrate_var.standard_name = "mole_concentration_of_nitrate_in_water_body"
nitrate_var.long_name = 'Nitrate concentration in water body'
nitrate_var.cell_methods = "time: mean"
nitrate_var.missing_value = -9999
nitrate_var.comment = "Reported in micromoles per liter (µmol/L). Equivalent to 1e-3 mol/m³, consistent with CF convention for mole concentrations. Measured in surface runoff from Rambla del Albujón, a freshwater input to the Mar Menor"
nitrate_var[:] = valores_con_nan

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
dataset.close()
# %%
ncfile = Dataset(f'{path}{nombre_fichero}.nc', "r")
txt_filename = f"{path}{nombre_fichero}.txt"

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
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')
"""
>=======>     >===>        >=>>=>   >=======>       >>       >===>>=====>     >===>
>=>         >=>    >=>   >=>    >=> >=>            >>=>           >=>       >=>    >=>
>=>       >=>        >=>  >=>       >=>           >> >=>          >=>     >=>        >=>
>=====>   >=>        >=>    >=>     >=====>      >=>  >=>         >=>     >=>        >=>
>=>       >=>        >=>       >=>  >=>         >=====>>=>        >=>     >=>        >=>
>=>         >=>     >=>  >=>    >=> >=>        >=>      >=>       >=>       >=>     >=>
>=>           >===>        >=>>=>   >=>       >=>        >=>      >=>         >===>

"""
nombre_fichero = 'CHS_ALBUJON_V3_PO4'
path= 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Runoff/CHS_ALBUJON_V3/'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title=f'{nombre_fichero}'
ncfile.institution="Confederación Hidrográfica del Segura (CHS)"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'CHS_ALBUJON_V3'
ncfile.project = 'Not associated with a specific project'
ncfile.source = 'In situ data collection'
ncfile.Conventions = "CF-1.8"

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
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values  # Se asigna directamente

valores_con_nan = (valor_array[:,2] *1000) / 94.97 #paso de mg L-1 a umol L-1 
valores_con_nan[np.isnan(valores_con_nan)] = -9999

fosfate_var = ncfile.createVariable('phosphate', np.float64, ('time',))
fosfate_var[:] = valores_con_nan
fosfate_var.standard_name = 'mole_concentration_of_phosphate_in_water_body'
fosfate_var.long_name ='Phosphate concentration in water body'
fosfate_var.cell_methods = "time: mean"
fosfate_var.missing_value = -9999
fosfate_var.comment = "Reported in micromoles per liter (µmol/L). Equivalent to 1e-3 mol/m³, consistent with CF convention for mole concentrations. Measured in surface runoff from Rambla del Albujón, a freshwater input to the Mar Menor"

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


# %%
dataset.close()
# %%
ncfile = Dataset(f'{path}{nombre_fichero}.nc', "r")
txt_filename = f"{path}{nombre_fichero}.txt"

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
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')

# %%