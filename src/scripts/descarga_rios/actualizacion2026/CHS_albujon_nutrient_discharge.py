from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import pandas as pd
import shutil
import sys
import os
import matplotlib.pyplot as plt
sys.path.append("../../")
from generar_txt import generar_txt
# %%
data = pd.read_csv('C:/Users/Julia/Documents/VSCODE_BELLICH/src/scripts/descarga_rios/actualizacion2026/datos/CHS/CHS_Albujon_water_chemistry.csv', sep=',', dtype={  
"Variable": "object",
"Unit":"object",
"Value":'float64', }, parse_dates= ['Date'],dayfirst=True )

for column in data.columns:
    print((data[column].iloc[0]))
    print(type(data[column].iloc[0]))
    print('-------------------------------')
# %%
#  conversión de fechas
# Fecha de referencia (1 de enero de 1970)

df_cond = data.loc[data['Variable'] == 'Electrical_Conductivity_20degC']
df_cond_media = (
    df_cond
    .groupby('Date', as_index=False)
    .agg({'Value': 'mean'})
)

print(df_cond_media.head())
df_cond = df_cond_media.copy()
# %%
# ANALISIS DUPLICADOS
# Número total de filas
print("Filas totales:", len(df_cond))
# Número de fechas únicas
print("Fechas únicas:", df_cond['Date'].nunique())
duplicadas = df_cond[df_cond.duplicated(subset='Date', keep=False)]
print(duplicadas.sort_values('Date'))
# -----------------------------------------------------------

fechas_unicas = np.sort(df_cond["Date"].unique())
fechas_unicas_ts = pd.to_datetime(fechas_unicas)

epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (fechas_unicas_ts - epoch) / pd.Timedelta(days=1)

# %%
nombre_fichero='CHS_Albujon_COND'
path= 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Runoff/CHS_Albujon/'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title=f'{nombre_fichero}'
ncfile.institution="Confederación Hidrográfica del Segura (CHS)"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'CHS'
ncfile.project = 'Not associated with a specific project'
ncfile.source = 'In situ data collection'
ncfile.Conventions = "CF-1.8"

# Cuantos parameter hay
lista_params = (data["Variable"].unique())
# 
ncfile.createDimension('time', len(dias_desde_1970))

for dim in ncfile.dimensions.items():
    print(dim)
# %%

time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:00"
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values  # Se asigna directamente

lat_var = ncfile.createVariable('latitude', np.float64, )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var.grid_mapping = "crs"
lat_var[:] =  37.716221

lon_var = ncfile.createVariable('longitude', np.float64, )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] = -0.861232 

# Lo paso a mS cm --> 1 mS/cm= 1000 μS/cm
valores_con_nan = df_cond['Value']
valores_con_nan[np.isnan(valores_con_nan)] = -9999 
# %%
conductivity_var = ncfile.createVariable('conductivity', np.float64, ('time',))
conductivity_var.units= 'uS cm-1' # unidades anteriores, no son estandar
# conductivity_var.units = 'mS cm-1' dejamos las unidades originales
conductivity_var.standard_name = "water_body_electrical_conductivity"
conductivity_var.long_name = 'Water body electrical conductivity'
conductivity_var.cell_methods= "time: mean" # no se si es mean o puntual
conductivity_var.missing_value = -9999 
conductivity_var.comment = 'Conductivity is measured in situ and corrected as equivalent at 20º. Conductivity probe model is not reported.' \
'Sampling hour is not recorded but usually it is the first light hours of the day'\
'Sampling intervals of variable length. From 2006 to 2016 an average of 6 sampling dates per year (2-13); from 2017 to  2021 sampling intervals were approximately fortnightly, and from 2022 onwards samplings are approximately weekly.'
conductivity_var[:] = valores_con_nan

crs = ncfile.createVariable("crs", "i")  # Dummy variable for coordinate reference system
crs.grid_mapping_name = "latitude_longitude"
crs.projection = "Geodetic"
crs.long_name = "WGS 84 / Geographic coordinates (EPSG:4326)"
crs.epsg_code = "EPSG:4326"
crs.semi_major_axis = 6378137.0
crs.inverse_flattening = 298.257223563
crs.comment = "Geographic coordinates are referenced to WGS 84 (EPSG:4326) in decimal degrees."

# %%
ncfile.close()

generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')

# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Runoff/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')

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
conductivity = dataset.variables["conductivity"][:]    #
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
fechas = np.array(fechas) 

# # %%
# Graficar cada parámetro en un subplot diferente
fig, ax = plt.subplots(1, 1, figsize=(10, 5 * (3)), sharex=True)
ax.plot(fechas , 
            conductivity,
            marker="o",
            linestyle="-")
    
ax.set_ylabel(f"conductividad ")
ax.legend()
ax.grid()

plt.xlabel("Fecha")
plt.tight_layout()
file_name = f'conductividad_temporal_profile.png'
plt.savefig(os.path.join(path, file_name), dpi=300)
plt.show()
# %%
dataset.close()
# %%


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
df_cond = data.loc[data['Variable'] == 'Nitrate']
df_cond_media = (
    df_cond
    .groupby('Date', as_index=False)
    .agg({'Value': 'mean'})
)

print(df_cond_media.head())
df_cond = df_cond_media.copy()
# %%
# ANALISIS DUPLICADOS
# Número total de filas
print("Filas totales:", len(df_cond))
# Número de fechas únicas
print("Fechas únicas:", df_cond['Date'].nunique())
duplicadas = df_cond[df_cond.duplicated(subset='Date', keep=False)]
print(duplicadas.sort_values('Date'))
# -----------------------------------------------------------

fechas_unicas = np.sort(df_cond["Date"].unique())
fechas_unicas_ts = pd.to_datetime(fechas_unicas)

epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (fechas_unicas_ts - epoch) / pd.Timedelta(days=1)

nombre_fichero = 'CHS_Albujon_NO3'
path= 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Runoff/CHS_Albujon/'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title=f'{nombre_fichero}'
ncfile.institution="Confederación Hidrográfica del Segura (CHS)"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'CHS'
ncfile.project = 'Not associated with a specific project'
ncfile.source = 'In situ data collection'
ncfile.Conventions = "CF-1.8"

# Cuantos parameter hay
lista_params = (data["Variable"].unique())
# 
ncfile.createDimension('time', len(dias_desde_1970))

for dim in ncfile.dimensions.items():
    print(dim)
# %%

time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:00"
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values  # Se asigna directamente

lat_var = ncfile.createVariable('latitude', np.float64, )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var.grid_mapping = "crs"
lat_var[:] =  37.716221

lon_var = ncfile.createVariable('longitude', np.float64, )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] = -0.861232 

# Lo paso a mS cm --> 1 mS/cm= 1000 μS/cm
valores_con_nan = df_cond['Value']
valores_con_nan[np.isnan(valores_con_nan)] = -9999 
# %%
conductivity_var = ncfile.createVariable('nitrate', np.float64, ('time',))
conductivity_var.units= 'mg L-1' # unidades qye dan en el excel
conductivity_var.standard_name = 'mass_concentration_of_nitrate_in_sea_water'
conductivity_var.long_name = 'Dissolved nitrate concentration in sea water'
conductivity_var.comment = 'Nitrate measured in the lab. Laboratory methods are not reported.'\
'Sampling hour is not recorded but usually it is the first light hours of the day'\
'Sampling intervals of variable length. From 2006 to 2016 an average of 6 sampling dates per year (2-13); from 2017 to 2021 sampling intervals were approximately fortnightly, and from 2022 onwards samplings are approximately weekly.'
conductivity_var[:] = valores_con_nan

crs = ncfile.createVariable("crs", "i")  # Dummy variable for coordinate reference system
crs.grid_mapping_name = "latitude_longitude"
crs.projection = "Geodetic"
crs.long_name = "WGS 84 / Geographic coordinates (EPSG:4326)"
crs.epsg_code = "EPSG:4326"
crs.semi_major_axis = 6378137.0
crs.inverse_flattening = 298.257223563
crs.comment = "Geographic coordinates are referenced to WGS 84 (EPSG:4326) in decimal degrees."

# %%
ncfile.close()

generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')

# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Runoff/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')

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
nitrate = dataset.variables["nitrate"][:]    # 
# fofate = dataset.variables["fosfate"][:]    #  
# nit = dataset.variables["conductivity"][:]    #
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
fechas = np.array(fechas) 

# # %%
# Graficar cada parámetro en un subplot diferente
fig, ax = plt.subplots(1, 1, figsize=(10, 5 * (3)), sharex=True)
ax.plot(fechas , 
            nitrate,
            marker="o",
            linestyle="-")
    
ax.set_ylabel(f"nitrate ")
ax.legend()
ax.grid()

plt.xlabel("Fecha")
plt.tight_layout()
file_name = f'nitrate_temporal_profile.png'
plt.savefig(os.path.join(path, file_name), dpi=300)
plt.show()
# %%
dataset.close()
"""
>=======>     >===>        >=>>=>   >=======>       >>       >===>>=====>     >===>
>=>         >=>    >=>   >=>    >=> >=>            >>=>           >=>       >=>    >=>
>=>       >=>        >=>  >=>       >=>           >> >=>          >=>     >=>        >=>
>=====>   >=>        >=>    >=>     >=====>      >=>  >=>         >=>     >=>        >=>
>=>       >=>        >=>       >=>  >=>         >=====>>=>        >=>     >=>        >=>
>=>         >=>     >=>  >=>    >=> >=>        >=>      >=>       >=>       >=>     >=>
>=>           >===>        >=>>=>   >=>       >=>        >=>      >=>         >===>

"""
df_cond = data.loc[data['Variable'] == 'Phosphate']
df_cond_media = (
    df_cond
    .groupby('Date', as_index=False)
    .agg({'Value': 'mean'})
)

print(df_cond_media.head())
df_cond = df_cond_media.copy()
# %%
# ANALISIS DUPLICADOS
# Número total de filas
print("Filas totales:", len(df_cond))
# Número de fechas únicas
print("Fechas únicas:", df_cond['Date'].nunique())
duplicadas = df_cond[df_cond.duplicated(subset='Date', keep=False)]
print(duplicadas.sort_values('Date'))
# -----------------------------------------------------------

fechas_unicas = np.sort(df_cond["Date"].unique())
fechas_unicas_ts = pd.to_datetime(fechas_unicas)

epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (fechas_unicas_ts - epoch) / pd.Timedelta(days=1)

nombre_fichero = 'CHS_Albujon_PO4'
path= 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Runoff/CHS_Albujon/'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title=f'{nombre_fichero}'
ncfile.institution="Confederación Hidrográfica del Segura (CHS)"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'CHS'
ncfile.project = 'Not associated with a specific project'
ncfile.source = 'In situ data collection'
ncfile.Conventions = "CF-1.8"

# Cuantos parameter hay
lista_params = (data["Variable"].unique())
# 
ncfile.createDimension('time', len(dias_desde_1970))

for dim in ncfile.dimensions.items():
    print(dim)
# %%

time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:00"
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values  # Se asigna directamente

lat_var = ncfile.createVariable('latitude', np.float64, )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var.grid_mapping = "crs"
lat_var[:] =  37.716221

lon_var = ncfile.createVariable('longitude', np.float64, )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] = -0.861232 

# Lo paso a mS cm --> 1 mS/cm= 1000 μS/cm
valores_con_nan = df_cond['Value']
valores_con_nan[np.isnan(valores_con_nan)] = -9999 
# %%
conductivity_var = ncfile.createVariable('phosphate', np.float64, ('time',))
conductivity_var.units= 'mg L-1' # unidades qye dan en el excel
conductivity_var.standard_name = 'mass_concentration_of_phosphate_in_sea_water'
conductivity_var.long_name = 'Dissolved phosphate concentration in sea water'
conductivity_var.comment = 'Phosphate measured in the lab. Laboratory methods are not reported.'\
'Sampling hour is not recorded but usually it is the first light hours of the day'\
'Sampling intervals of variable length. From 2006 to 2016 an average of 6 sampling dates per year (2-13); from 2017 to 2021 sampling intervals were approximately fortnightly, and from 2022 onwards samplings are approximately weekly.'
conductivity_var[:] = valores_con_nan

crs = ncfile.createVariable("crs", "i")  # Dummy variable for coordinate reference system
crs.grid_mapping_name = "latitude_longitude"
crs.projection = "Geodetic"
crs.long_name = "WGS 84 / Geographic coordinates (EPSG:4326)"
crs.epsg_code = "EPSG:4326"
crs.semi_major_axis = 6378137.0
crs.inverse_flattening = 298.257223563
crs.comment = "Geographic coordinates are referenced to WGS 84 (EPSG:4326) in decimal degrees."

# %%
ncfile.close()

generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')

# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Runoff/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')

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
fofate = dataset.variables["phosphate"][:]    #  
# nit = dataset.variables["conductivity"][:]    #
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
fechas = np.array(fechas) 

# # %%
# Graficar cada parámetro en un subplot diferente
fig, ax = plt.subplots(1, 1, figsize=(10, 5 * (3)), sharex=True)
ax.plot(fechas , 
            fofate,
            marker="o",
            linestyle="-")
    
ax.set_ylabel(f"fosfate ")
ax.legend()
ax.grid()

plt.xlabel("Fecha")
plt.tight_layout()
file_name = f'fosfato_temporal_profile.png'
plt.savefig(os.path.join(path, file_name), dpi=300)
plt.show()
# %%
dataset.close()
# %%