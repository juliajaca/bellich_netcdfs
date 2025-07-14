from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
from generar_txt import generar_txt
# %%
data = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/rios/Albujon_discharge_daily.xlsx', dtype={  
"Flow_mean": "float64", "Discharge":"float64" }, parse_dates= ['Date'], sheet_name='Data_CHS')

# %%
#  conversión de fechas : Fecha de referencia (1 de enero de 1970)
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (data.Date - epoch) / pd.Timedelta(days=1)

# %%
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/AportesContinentales/CHS_ALBUJON_v1/CHS_ALBUJON_V1_Q_MEAN_DAILY.nc'

ncfile = Dataset(path, mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='Albujon flow daily'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'
ncfile.source = 'Confederación Hidrográfica del Segura (CHS)'
ncfile.comment = (
    "Flow measured every 5 minutes by CHS automatic gauge; daily mean of flows.\n"
    "Stream gauge located at the river mouth."
)
ncfile.Conventions = "CF-1.8"

# Crear dimensiones
ncfile.createDimension('time', len(data))

for dim in ncfile.dimensions.items():
    print(dim)
# %%
date_var = ncfile.createVariable('time', np.float64, ('time'))
date_var.units= "days since 1970-01-01 00:00:0"
date_var.standard_name = "time"
date_var.calendar = 'gregorian'
date_var[:] = dias_desde_1970.values

flow_var = ncfile.createVariable('water_volume_transport', np.float64, ('time',))
flow_var.units = 'm3 s-1'
flow_var.standard_name = "water_volume_transport_in_river_channel"
flow_var.long_name = "Mean river discharge"
flow_var[:] = data["Flow_mean"].values

lat_var = ncfile.createVariable('latitude', np.float64, )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var[:] = 37.716068

lon_var = ncfile.createVariable('longitude', np.float64, )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var[:] =   -0.861019 

crs = ncfile.createVariable("crs", "i")  # Dummy variable for coordinate reference system
crs.long_name=  "Reference coordinate system"
crs.grid_mapping_name = "latitude_longitude"
crs.epsg_code = "EPSG:4326"

# %%
ncfile.close()

# %%
txt_filename = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/AportesContinentales/CHS_ALBUJON_v1/CHS_ALBUJON_V1_Q_MEAN_DAILY.txt'
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
# DISCHARGE
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/AportesContinentales/CHS_ALBUJON_v1/CHS_ALBUJON_V1_Q_ACCUM_DAILY.nc'

ncfile = Dataset(path, mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='Albujon discharge daily'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'
ncfile.source = 'Confederación Hidrográfica del Segura (CHS)'
ncfile.comment = (
    "The discharges are calculated by integrating the instantaneous flow over the day.\n"
    "Stream gauge located at the river mouth."
)
ncfile.Conventions = "CF-1.8"

# Crear dimensiones
ncfile.createDimension('time', len(data))

for dim in ncfile.dimensions.items():
    print(dim)
# %%
date_var = ncfile.createVariable('time', np.float64, ('time'))
date_var.units= "days since 1970-01-01 00:00:0"
date_var.standard_name = "time"
date_var.calendar = 'gregorian'
date_var[:] = dias_desde_1970.values

disch_var = ncfile.createVariable('accumulated_water_volume_transport', np.float64, ('time',))
disch_var.units = 'm3'
disch_var.long_name = "Accumulated daily discharge"
disch_var[:] = data["Discharge"].values

lat_var = ncfile.createVariable('latitude', np.float64, )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var[:] = 37.716068

lon_var = ncfile.createVariable('longitude', np.float64, )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var[:] =   -0.861019 

crs = ncfile.createVariable("crs", "i")  # Dummy variable for coordinate reference system
crs.long_name=  "Reference coordinate system"
crs.grid_mapping_name = "latitude_longitude"
crs.epsg_code = "EPSG:4326"

# %%
ncfile.close()

# %%

txt_filename = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/AportesContinentales/CHS_ALBUJON_v1/CHS_ALBUJON_V1_Q_ACCUM_DAILY.txt'
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
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/AportesContinentales/CARM_ALBUJON/CARM_ALBUJON_Q_MEAN_DAILY.nc'

ncfile = Dataset(path, mode='w', format='NETCDF3_CLASSIC')
print(ncfile)
# # %% CREAR ATRIBUTOS GLOBALES
ncfile.title='Albujon flow daily'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'
ncfile.source = 'https://canalmarmenor.carm.es/monitorizacion/monitorizacion-de-parametros/aforos/'
ncfile.comment = (
    "Flow measured at a single point in time on that date.\n"
    "Stream gauge located at the river mouth."
)
ncfile.Conventions = "CF-1.8"

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
lat_var[:] = 37.716068

lon_var = ncfile.createVariable('longitude', np.float64, )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var[:] =   -0.861019 

crs = ncfile.createVariable("crs", "i")  # Dummy variable for coordinate reference system
crs.long_name=  "Reference coordinate system"
crs.grid_mapping_name = "latitude_longitude"
crs.epsg_code = "EPSG:4326"

flow_var2 = ncfile.createVariable('water_volume_transport', np.float64, ('time',))
flow_var2.units = 'm3 s-1'
flow_var2.standard_name = "water_volume_transport_in_river_channel"
flow_var2.long_name = "Mean river discharge"
flow_var2[:] = data["Flow"].values

ncfile.close()


# %%

txt_filename = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/AportesContinentales/CARM_ALBUJON/CARM_ALBUJON_Q_MEAN_DAILY.txt'
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
# %% CARM 2
# %%
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/AportesContinentales/CARM_ALBUJON/CARM_ALBUJON_Q_ACCUM_DAILY.nc'

ncfile = Dataset(path, mode='w', format='NETCDF3_CLASSIC')
print(ncfile)
# # %% CREAR ATRIBUTOS GLOBALES
ncfile.title='Albujon discharge daily'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'
ncfile.source = 'https://canalmarmenor.carm.es/monitorizacion/monitorizacion-de-parametros/aforos/'
ncfile.comment = (
    "The discharges are calculated by integrating the instantaneous flow over the day.\n"
    "Stream gauge located at the river mouth."
)
ncfile.Conventions = "CF-1.8"

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
lat_var[:] = 37.716068

lon_var = ncfile.createVariable('longitude', np.float64, )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var[:] =   -0.861019 

crs = ncfile.createVariable("crs", "i")  # Dummy variable for coordinate reference system
crs.long_name=  "Reference coordinate system"
crs.grid_mapping_name = "latitude_longitude"
crs.epsg_code = "EPSG:4326"

disch_var2 = ncfile.createVariable('accumulated_water_volume_transport', np.float64, ('time',))
disch_var2.units = 'm³'
disch_var2.long_name = "Accumulated daily discharge"
disch_var2[:] = data["Discharge"].values

ncfile.close()

# %%
txt_filename = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/AportesContinentales/CARM_ALBUJON/CARM_ALBUJON_Q_ACCUM_DAILY.txt'
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