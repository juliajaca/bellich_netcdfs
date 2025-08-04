from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
from generar_txt import generar_txt
# %%
data = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/rios/Albujon_discharge_daily.xlsx', dtype={  
"Flow_mean": "float64", "Discharge":"float64" }, parse_dates= ['Date'], sheet_name='Data_CHS')

nombre_fichero= 'CHS_ALBUJON_V1_Q_MEAN_DAILY'
# %%
#  conversión de fechas : Fecha de referencia (1 de enero de 1970)
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (data.Date - epoch) / pd.Timedelta(days=1)

# %%
path = f'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Runoff/CHS_ALBUJON_V1/{nombre_fichero}.nc'

ncfile = Dataset(path, mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title=f'{nombre_fichero}'
ncfile.institution="Confederación Hidrográfica del Segura (CHS)"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'CHS_ALBUJON_V1'
ncfile.project = 'Not associated with a specific project'
ncfile.source = 'In situ data collection'
ncfile.Conventions = "CF-1.8"
ncfile.comment = (
    "Flow measured every 5 minutes by CHS automatic gauge; daily mean of flows.\n"
    "Stream gauge located at the river mouth."
)

# Crear dimensiones
ncfile.createDimension('time', len(data))

for dim in ncfile.dimensions.items():
    print(dim)
# %%
date_var = ncfile.createVariable('time', np.float64, ('time'))
date_var.units= "days since 1970-01-01 00:00:0"
date_var.calendar = 'gregorian'
date_var.standard_name = "time"
date_var[:] = dias_desde_1970.values

lat_var = ncfile.createVariable('latitude', np.float64, )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var.grid_mapping = "crs"
lat_var[:] = 37.716068

lon_var = ncfile.createVariable('longitude', np.float64, )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] =   -0.861019 

flow_var = ncfile.createVariable('water_volume_transport', np.float64, ('time',))
flow_var.units = 'm3 s-1'
flow_var.standard_name = "volume_transport_in_river_channel"
flow_var.long_name = 'Daily Mean River Discharge'
valores_con_nan = data["Flow_mean"].values
valores_con_nan[np.isnan(valores_con_nan)] = -9999 
flow_var.cell_methods= "time: mean"
flow_var[:] =  valores_con_nan
flow_var.missing_value = -9999 
flow_var.grid_mapping = "crs"

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


# %%
txt_filename = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Runoff/CHS_ALBUJON_V1/CHS_ALBUJON_V1_Q_MEAN_DAILY.txt'
ncfile = Dataset(path, "r")

with open(txt_filename, "w") as f:
    # Escribir el formato
    f.write(f"Format:\n\t{ncfile.file_format}\n")
    f.write("\nGlobal Attributes:\n")
    for attr in ncfile.ncattrs():
        f.write(f"\t{attr:<15} = '{getattr(ncfile, attr)}'\n")
    f.write("\nDimensions:\n")
    for dim_name, dim in ncfile.dimensions.items():
        f.write(f"\t{dim_name:<15} = {len(dim)}\n")
    f.write("\nVariables:\n")
    for var_name, var in ncfile.variables.items():
        # Obtener dimensiones y tipo de datos
        dimensions = ", ".join(var.dimensions)
        dtype = var.dtype
        f.write(f"    {var_name:<15}\n")
        f.write(f"\tSize: {var.shape}\n")
        f.write(f"\tDimensions: {dimensions}\n")
        f.write(f"\tDatatype: {dtype}\n")
        if var.ncattrs():
            f.write("\tAttributes:\n")
            for attr in var.ncattrs():
                f.write(f"\t\t{attr:<15} = '{getattr(var, attr)}'\n")
ncfile.close()

print(f"Archivo '{txt_filename}' generado con éxito.")

# %%
#  COPIAR
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Runoff/CHS_ALBUJON_V1/'
shutil.copy(path,f'{ruta_destino}{nombre_fichero}.nc')



# DISCHARGE
nombre_fichero= 'CHS_ALBUJON_V1_Q_ACCUM_DAILY'
path = f'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Runoff/CHS_ALBUJON_V1/{nombre_fichero}.nc'

ncfile = Dataset(path, mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title=f'{nombre_fichero}'
ncfile.institution="Confederación Hidrográfica del Segura (CHS)"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'CHS_ALBUJON_V1'
ncfile.project = 'Not associated with a specific project'
ncfile.source = 'In situ data collection'
ncfile.Conventions = "CF-1.8"
ncfile.comment = (
    "Stream gauge located at the river mouth.\n" \
    "The discharges are calculated integrating the instantaneous flow for the whole day"
)

# Crear dimensiones
ncfile.createDimension('time', len(data))

for dim in ncfile.dimensions.items():
    print(dim)
# %%
date_var = ncfile.createVariable('time', np.float64, ('time'))
date_var.units= "days since 1970-01-01 00:00:0"
date_var.calendar = 'gregorian'
date_var.standard_name = "time"
date_var[:] = dias_desde_1970.values

lat_var = ncfile.createVariable('latitude', np.float64, )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var.grid_mapping = "crs"
lat_var[:] = 37.716068

lon_var = ncfile.createVariable('longitude', np.float64, )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] =   -0.861019 

valores_con_nan= data["Discharge"].values
valores_con_nan[np.isnan(valores_con_nan)] = -9999 

disch_var = ncfile.createVariable('accumulated_water_volume_transport', np.float64, ('time',))
disch_var.units = 'm3'
disch_var.standard_name =' volume_transport_in_river_channel'
disch_var.long_name= 'Daily Accumulated River Discharge'
disch_var.cell_methods= 'time: sum'
disch_var.missing_value = -9999 
disch_var.grid_mapping = "crs"
disch_var[:] = valores_con_nan

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

# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Runoff/CHS_ALBUJON_V1/'
shutil.copy(path,f'{ruta_destino}{nombre_fichero}.nc')

txt_filename = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Runoff/CHS_ALBUJON_v1/CHS_ALBUJON_V1_Q_ACCUM_DAILY.txt'
ncfile = Dataset(path, "r")

with open(txt_filename, "w") as f:
    # Escribir el formato
    f.write(f"Format:\n\t{ncfile.file_format}\n")
    f.write("\nGlobal Attributes:\n")
    for attr in ncfile.ncattrs():
        f.write(f"\t{attr:<15} = '{getattr(ncfile, attr)}'\n")
    f.write("\nDimensions:\n")
    for dim_name, dim in ncfile.dimensions.items():
        f.write(f"\t{dim_name:<15} = {len(dim)}\n")
    f.write("\nVariables:\n")
    for var_name, var in ncfile.variables.items():
        # Obtener dimensiones y tipo de datos
        dimensions = ", ".join(var.dimensions)
        dtype = var.dtype
        f.write(f"    {var_name:<15}\n")
        f.write(f"\tSize: {var.shape}\n")
        f.write(f"\tDimensions: {dimensions}\n")
        f.write(f"\tDatatype: {dtype}\n")
        if var.ncattrs():
            f.write("\tAttributes:\n")
            for attr in var.ncattrs():
                f.write(f"\t\t{attr:<15} = '{getattr(var, attr)}'\n")
ncfile.close()

print(f"Archivo '{txt_filename}' generado con éxito.")

# %%
# %%
"""
    >=>          >>       >======>     >=>       >=>
 >=>   >=>      >>=>      >=>    >=>   >> >=>   >>=>
>=>            >> >=>     >=>    >=>   >=> >=> > >=>
>=>           >=>  >=>    >> >==>      >=>  >=>  >=>
>=>          >=====>>=>   >=>  >=>     >=>   >>  >=>
 >=>   >=>  >=>      >=>  >=>    >=>   >=>       >=>
   >===>   >=>        >=> >=>      >=> >=>       >=>

"""
data = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/rios/Albujon_discharge_daily.xlsx', dtype={  
"Flow":'float64', 'Discharge':'float64'}, parse_dates= ['Date'], sheet_name='Data_CARM')
# %%
#  conversión de fechas Fecha de referencia (1 de enero de 1970)
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (data.Date - epoch) / pd.Timedelta(days=1)

# %%
nombre_fichero= 'CARM_ALBUJON_Q_DAILY'
path = f'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Runoff/CARM_ALBUJON/'

ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')

print(ncfile)
# # %% CREAR ATRIBUTOS GLOBALES
ncfile.title = nombre_fichero 
ncfile.institution="CARM"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'CHS_ALBUJON_V1'
ncfile.project = 'Not associated with a specific project'
ncfile.source = 'In situ data collection https://canalmarmenor.carm.es/monitorizacion/monitorizacion-de-parametros/aforos/'
ncfile.Conventions = "CF-1.8"
ncfile.comment = (
    "Flow measured in a time point in time on that date.\n"
    "Stream gauge located at the river mouth."
)

# %%
ncfile.createDimension('time', len(data))
date_var = ncfile.createVariable('time', np.float64, ('time'))
date_var.units= "days since 1970-01-01 00:00:0"
date_var.standard_name = "time"
date_var.calendar = 'gregorian'
date_var[:] = dias_desde_1970.values

lat_var = ncfile.createVariable('latitude', np.float64, )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var.grid_mapping = "crs"
lat_var[:] = 37.716068

lon_var = ncfile.createVariable('longitude', np.float64, )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] = -0.861019 

valores_con_nan = data["Flow"].values
valores_con_nan[np.isnan(valores_con_nan)] = -9999 
flow_var2 = ncfile.createVariable('water_volume_transport', np.float64, ('time',))
flow_var2.units = 'm3 s-1'
flow_var2.standard_name = "volume_transport_in_river_channel"
flow_var2.long_name = 'Data Point River Discharge'
flow_var2.cell_methods= 'time: point'
flow_var2.missing_value = -9999 
flow_var2.grid_mapping = "crs"
flow_var2[:] = valores_con_nan

crs = ncfile.createVariable("crs", "i")  # Dummy variable for coordinate reference system
crs.grid_mapping_name = "latitude_longitude"
crs.projection = "Geodetic"
crs.long_name = "WGS 84 / Geographic coordinates (EPSG:4326)"
crs.epsg_code = "EPSG:4326"
crs.semi_major_axis = 6378137.0
crs.inverse_flattening = 298.257223563
crs.comment = "Geographic coordinates are referenced to WGS 84 (EPSG:4326) in decimal degrees."

ncfile.close()
# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Runoff/CARM_ALBUJON/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')

# %%

txt_filename = f'{path}{nombre_fichero}.txt'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', "r")

with open(txt_filename, "w") as f:
    # Escribir el formato
    f.write(f"Format:\n\t{ncfile.file_format}\n")
    f.write("\nGlobal Attributes:\n")
    for attr in ncfile.ncattrs():
        f.write(f"\t{attr:<15} = '{getattr(ncfile, attr)}'\n")
    f.write("\nDimensions:\n")
    for dim_name, dim in ncfile.dimensions.items():
        f.write(f"\t{dim_name:<15} = {len(dim)}\n")
    f.write("\nVariables:\n")
    for var_name, var in ncfile.variables.items():
        # Obtener dimensiones y tipo de datos
        dimensions = ", ".join(var.dimensions)
        dtype = var.dtype
        f.write(f"    {var_name:<15}\n")
        f.write(f"\tSize: {var.shape}\n")
        f.write(f"\tDimensions: {dimensions}\n")
        f.write(f"\tDatatype: {dtype}\n")
        if var.ncattrs():
            f.write("\tAttributes:\n")
            for attr in var.ncattrs():
                f.write(f"\t\t{attr:<15} = '{getattr(var, attr)}'\n")
ncfile.close()

print(f"Archivo '{txt_filename}' generado con éxito.")
# %% CARM 2
# %%
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Runoff/CARM_ALBUJON/'
nombre_fichero= 'CARM_ALBUJON_Q_ACCUM_DAILY'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)
# # %% CREAR ATRIBUTOS GLOBALES
ncfile.title = nombre_fichero 
ncfile.institution="CARM"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'CHS_ALBUJON_V1'
ncfile.project = 'Not associated with a specific project'
ncfile.source = 'In situ data collection https://canalmarmenor.carm.es/monitorizacion/monitorizacion-de-parametros/aforos/'
ncfile.Conventions = "CF-1.8"
ncfile.comment = (
    "The discharges are calculated by integrating the instantaneous flow over the day.\n"
    "Stream gauge located at the river mouth."
)

ncfile.createDimension('time', len(data))

# %%
date_var = ncfile.createVariable('time', np.float64, ('time'))
date_var.units= "days since 1970-01-01 00:00:0"
date_var.standard_name = "time"
date_var.calendar = 'gregorian'
date_var[:] = dias_desde_1970.values

lat_var = ncfile.createVariable('latitude', np.float64, )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var.grid_mapping = "crs"
lat_var[:] = 37.716068

lon_var = ncfile.createVariable('longitude', np.float64, )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] =   -0.861019 

disch_var2 = ncfile.createVariable('accumulated_water_volume_transport', np.float64, ('time',))
disch_var2.units = 'm³'
disch_var2.standard_name =' volume_transport_in_river_channel'
disch_var2.long_name= 'Daily Accumulated River Discharge'
disch_var2.cell_methods= 'time: sum'
disch_var2.missing_value = -9999 
disch_var2.grid_mapping = "crs"

valores_con_nan= data["Discharge"].values
valores_con_nan[np.isnan(valores_con_nan)] = -9999 
disch_var2[:] = valores_con_nan

crs = ncfile.createVariable("crs", "i")  # Dummy variable for coordinate reference system
crs.grid_mapping_name = "latitude_longitude"
crs.projection = "Geodetic"
crs.long_name = "WGS 84 / Geographic coordinates (EPSG:4326)"
crs.epsg_code = "EPSG:4326"
crs.semi_major_axis = 6378137.0
crs.inverse_flattening = 298.257223563
crs.comment = "Geographic coordinates are referenced to WGS 84 (EPSG:4326) in decimal degrees."

ncfile.close()
# %%
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')
# %%
txt_filename = f'{path}{nombre_fichero}.txt'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', "r")

with open(txt_filename, "w") as f:
    # Escribir el formato
    f.write(f"Format:\n\t{ncfile.file_format}\n")
    f.write("\nGlobal Attributes:\n")
    for attr in ncfile.ncattrs():
        f.write(f"\t{attr:<15} = '{getattr(ncfile, attr)}'\n")
    f.write("\nDimensions:\n")
    for dim_name, dim in ncfile.dimensions.items():
        f.write(f"\t{dim_name:<15} = {len(dim)}\n")
    f.write("\nVariables:\n")
    for var_name, var in ncfile.variables.items():
        # Obtener dimensiones y tipo de datos
        dimensions = ", ".join(var.dimensions)
        dtype = var.dtype
        f.write(f"    {var_name:<15}\n")
        f.write(f"\tSize: {var.shape}\n")
        f.write(f"\tDimensions: {dimensions}\n")
        f.write(f"\tDatatype: {dtype}\n")
        if var.ncattrs():
            f.write("\tAttributes:\n")
            for attr in var.ncattrs():
                f.write(f"\t\t{attr:<15} = '{getattr(var, attr)}'\n")
ncfile.close()

print(f"Archivo '{txt_filename}' generado con éxito.")


# %%
# %%
# Hacer los plots
# plt.figure(figsize=(10, 5))
# plt.plot(fechas, flujo, marker="o", linestyle="-", color="b", label="Flujo (m³/s)")
# plt.xlabel("Fecha"); plt.ylabel("Flujo (m³/s)"); plt.title("Serie Temporal del FLOW CHS"); plt.xticks(rotation=45); plt.legend(); plt.grid(); plt.show()
# # %%
# plt.figure(figsize=(10, 5))
# plt.plot(fechas, disch, marker="o", linestyle="-", color="b", label="Disch (m³)")
# plt.xlabel("Fecha");plt.ylabel("Discharge (m³)");plt.title("Serie Temporal del Discharge CHS");plt.xticks(rotation=45);plt.legend();plt.grid();plt.show()
# # %%
# plt.figure(figsize=(10, 5))
# plt.plot(fechas, flujo2, marker="o", linestyle="-", color="b", label="FLOW (m³/s)")
# plt.xlabel("Fecha");plt.ylabel("FLOW (m³/s)");plt.title("Serie Temporal del Flow carm");plt.xticks(rotation=45);plt.legend();plt.grid();plt.show()
# # %%
# plt.figure(figsize=(10, 5))
# plt.plot(fechas, disch2, marker="o", linestyle="-", color="b", label="Discharge (m³)")
# plt.xlabel("Fecha");plt.ylabel("Discharge (m³)");plt.title("Serie Temporal del discharge carm");plt.xticks(rotation=45);plt.legend();plt.grid();plt.show()