from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# %%
data = pd.read_excel('../datos/rios/Albujon_discharge_daily.xlsx', dtype={  
"Flow_CHS": "float32", "Disch_CHS":"float32", "Flow_CARM":'float32', 'Disch_CARM':'float32', 'Dif_disch':'float32' }, parse_dates= ['Date'], sheet_name='Combined_Data')
# %%
for column in data.columns:
    print((data[column].iloc[0]))
    print(type(data[column].iloc[0]))
# %%
#  conversión de fechas
# Fecha de referencia (1 de enero de 1970)
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (data.Date - epoch) / pd.Timedelta(days=1)


"""
  **
 ***
//**
 /**
 /**
 /**
 ****
////
"""
# %%
ncfile = Dataset('albujon_discharge_daily_flow_chs.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='Albujon flow daily CHS'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.poject = 'BELLICH'
ncfile.source = 'Confederación Hidrográfica del Segura (CHS)'
ncfile.conventions = 'XXX'
ncfile.coordinates = 'ETRS89 / UTM zone 30N --> EPSG:25830: [688529.72,4176466.36] \n' \
'ETRS89 --> EPSG:4258: [-0.861019,37.716068] \n' \
'WGS84 --> EPSG:4326: [-0.861019,37.716068]'
ncfile.notes= "Flow measured each 5' by CHS authomatic gauge; daily mean of flows. Stream gauge in the mouth.  The discharges are calculated integrating the instantaneous flow for the whole day"
# Crear dimensiones
date_dim = ncfile.createDimension('date', len(data))
flow_dim = ncfile.createDimension('flow_chs', len(data))
for dim in ncfile.dimensions.items():
    print(dim)
# %%
date_var = ncfile.createVariable('date', np.float32, ('date'))
date_var.long_name = 'Tiempo en días desde 1 enero 1970'
date_var.units = 'day'
date_var[:] = dias_desde_1970.values

flow_var = ncfile.createVariable('flow_chs', np.float32, ('date',))
flow_var.units = 'm³/s'
flow_var.long_name = 'Flujo medio en m3/s'
flow_var[:] = data["Flow_CHS"].values
# %%
ncfile.close()
# %% COMPROBACION
dataset = Dataset('albujon_discharge_daily_flow_chs.nc', "r")
print(dataset)
print(dataset.variables.keys())  # Ver las variables en el archivo
print('los fechas data son')
discharge_data = dataset.variables["date"][:]
print(discharge_data)
print('los flows son ')
print(dataset.variables['flow_chs'][:])

# print("\n🔹 Atributos de las Variables:")
# for var_name in dataset.variables:
#     print(f"\nVariable: {var_name}")
#     for attr in dataset.variables[var_name].ncattrs():
#         print(f"  {attr}: {dataset.variables[var_name].getncattr(attr)}")

# print("\n🔹 Atributos Globales:")
# for attr in dataset.ncattrs():
#     print(f"{attr}: {dataset.getncattr(attr)}")

# %%
# Leer las variables
tiempo = dataset.variables["date"][:]  # Días desde 1970
flujo = dataset.variables["flow_chs"][:]    # Flujo en m3/s

# Convertir tiempo a formato datetime para mejor visualización
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()

# Hacer el plot
plt.figure(figsize=(10, 5))
plt.plot(fechas, flujo, marker="o", linestyle="-", color="b", label="Flujo (m³/s)")
plt.xlabel("Fecha"); plt.ylabel("Flujo (m³/s)"); plt.title("Serie Temporal del FLOW CHS"); plt.xticks(rotation=45); plt.legend(); plt.grid(); plt.show()
# %%
"""
  ****
 */// *
/    /*
   ***
  *//
 *
/******
//////
"""
# %%
ncfile = Dataset('albujon_discharge_daily_discharge_chs.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='Albujon discharge daily CHS'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.poject = 'BELLICH'
ncfile.source = 'Confederación Hidrográfica del Segura (CHS)'
ncfile.conventions = 'XXX'
ncfile.coordinates = 'ETRS89 / UTM zone 30N --> EPSG:25830: [688529.72,4176466.36] \n' \
'ETRS89 --> EPSG:4258: [-0.861019,37.716068] \n' \
'WGS84 --> EPSG:4326: [-0.861019,37.716068]'
ncfile.notes= "Flow measured each 5' by CHS authomatic gauge; daily mean of flows. Stream gauge in the mouth.  The discharges are calculated integrating the instantaneous flow for the whole day"
# Crear dimensiones
date_dim = ncfile.createDimension('date', len(data))
disch_dim = ncfile.createDimension('disch_chs', len(data))
for dim in ncfile.dimensions.items():
    print(dim)
# %%
date_var = ncfile.createVariable('date', np.float32, ('date'))
date_var.long_name = 'Tiempo en días desde 1 enero 1970'
date_var.units = 'day'
date_var[:] = dias_desde_1970.values

flow_var = ncfile.createVariable('disch_chs', np.float32, ('date',))
flow_var.units = 'm³'
flow_var.long_name = 'Accumulated daily discharge'
flow_var[:] = data["Disch_CHS"].values
# %%
ncfile.close()
# %% COMPROBACION
dataset = Dataset('albujon_discharge_daily_discharge_chs.nc', "r")
print(dataset)
print(dataset.variables.keys())  # Ver las variables en el archivo
print('los fechas data son')
print(dataset.variables["date"][:])
print('los disharge  son ')
print(dataset.variables['disch_chs'][:])

# %% Leer las variables
tiempo = dataset.variables["date"][:]  # Días desde 1970
flujo = dataset.variables["disch_chs"][:]    # Flujo en m3
# Convertir tiempo a formato datetime para mejor visualización
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()
# Hacer el plot
plt.figure(figsize=(10, 5))
plt.plot(fechas, flujo, marker="o", linestyle="-", color="b", label="Disch (m³)")
plt.xlabel("Fecha");plt.ylabel("Discharge (m³)");plt.title("Serie Temporal del Discharge CHS");plt.xticks(rotation=45);plt.legend();plt.grid();plt.show()
# %% 
"""
  ****
 */// *
/    /*
   ***
  /// *
 *   /*
/ ****
 ////
"""
# %%
ncfile = Dataset('albujon_discharge_daily_flow_carm.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='Albujon flow daily CARM'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.poject = 'BELLICH'
ncfile.source = 'https://canalmarmenor.carm.es/monitorizacion/monitorizacion-de-parametros/aforos/'
ncfile.conventions = 'XXX'
ncfile.coordinates = 'ETRS89 / UTM zone 30N --> EPSG:25830: [688529.72,4176466.36] \n' \
'ETRS89 --> EPSG:4258: [-0.861019,37.716068] \n' \
'WGS84 --> EPSG:4326: [-0.861019,37.716068]'
ncfile.notes= "Flow measured in a time point that date. The discharges are calculated integrating the instantaneous flow for the whole day"
# Crear dimensiones
ncfile.createDimension('date', len(data))
ncfile.createDimension('flow_carm', len(data))
for dim in ncfile.dimensions.items():
    print(dim)
# %%
date_var = ncfile.createVariable('date', np.float32, ('date'))
date_var.long_name = 'Tiempo en días desde 1 enero 1970'
date_var.units = 'day'
date_var[:] = dias_desde_1970.values

flow_var = ncfile.createVariable('flow_carm', np.float32, ('date',))
flow_var.units = 'm³/s'
flow_var.long_name = 'Mean flow'
flow_var[:] = data["Flow_CARM"].values
# %%
ncfile.close()
# %% COMPROBACION
dataset = Dataset('albujon_discharge_daily_flow_carm.nc', "r")
print(dataset)
print(dataset.variables.keys())  # Ver las variables en el archivo
print('los fechas data son')
print(dataset.variables["date"][:])
print('los flow  son ')
print(dataset.variables['flow_carm'][:])

# %% Leer las variables
tiempo = dataset.variables["date"][:]  # Días desde 1970
flujo = dataset.variables["flow_carm"][:]    # Flujo en m3
# Convertir tiempo a formato datetime para mejor visualización
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()
# Hacer el plot
plt.figure(figsize=(10, 5))
plt.plot(fechas, flujo, marker="o", linestyle="-", color="b", label="FLOW (m³/s)")
plt.xlabel("Fecha");plt.ylabel("FLOW (m³/s)");plt.title("Serie Temporal del Flow carm");plt.xticks(rotation=45);plt.legend();plt.grid();plt.show()

# %%
"""
    **
   */*
  * /*
 ******
/////*
    /*
    /*
    /
# """
# %%
ncfile = Dataset('albujon_discharge_daily_discharge_carm.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='Albujon discharge daily CARM'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.poject = 'BELLICH'
ncfile.source = 'https://canalmarmenor.carm.es/monitorizacion/monitorizacion-de-parametros/aforos/'
ncfile.conventions = 'XXX'
ncfile.coordinates = 'ETRS89 / UTM zone 30N --> EPSG:25830: [688529.72,4176466.36] \n' \
'ETRS89 --> EPSG:4258: [-0.861019,37.716068] \n' \
'WGS84 --> EPSG:4326: [-0.861019,37.716068]'
ncfile.notes= "Flow measured in a time point that date. The discharges are calculated integrating the instantaneous flow for the whole day"
# Crear dimensiones
ncfile.createDimension('date', len(data))
ncfile.createDimension('disch_carm', len(data))
for dim in ncfile.dimensions.items():
    print(dim)
# %%
date_var = ncfile.createVariable('date', np.float32, ('date'))
date_var.long_name = 'Tiempo en días desde 1 enero 1970'
date_var.units = 'day'
date_var[:] = dias_desde_1970.values

flow_var = ncfile.createVariable('disch_carm', np.float32, ('date',))
flow_var.units = 'm³'
flow_var.long_name = 'Accumulated daily discharge'
flow_var[:] = data["Disch_CARM"].values
# %%
ncfile.close()
# %% COMPROBACION
dataset = Dataset('albujon_discharge_daily_discharge_carm.nc', "r")
print(dataset)
print(dataset.variables.keys())  # Ver las variables en el archivo
print('los fechas data son')
print(dataset.variables["date"][:])
print('los flow  son ')
print(dataset.variables['disch_carm'][:])

# %% Leer las variables
tiempo = dataset.variables["date"][:]  # Días desde 1970
flujo = dataset.variables["disch_carm"][:]    # Flujo en m3
# Convertir tiempo a formato datetime para mejor visualización
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()
# Hacer el plot
plt.figure(figsize=(10, 5))
plt.plot(fechas, flujo, marker="o", linestyle="-", color="b", label="Discharge (m³)")
plt.xlabel("Fecha");plt.ylabel("Discharge (m³)");plt.title("Serie Temporal del discharge carm");plt.xticks(rotation=45);plt.legend();plt.grid();plt.show()
# %%
"""
 ******
/*////
/*****
///// *
     /*
 *   /*
/ ****
 ////
"""
# %%
ncfile = Dataset('albujon_discharge_daily_dif_disch.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='Albujon dif discharge daily'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.poject = 'BELLICH'
ncfile.source = 'https://canalmarmenor.carm.es/monitorizacion/monitorizacion-de-parametros/aforos/ y Confederación Hidrográfica del Segura (CHS)'
ncfile.conventions = 'XXX'
ncfile.coordinates = 'ETRS89 / UTM zone 30N --> EPSG:25830: [688529.72,4176466.36] \n' \
'ETRS89 --> EPSG:4258: [-0.861019,37.716068] \n' \
'WGS84 --> EPSG:4326: [-0.861019,37.716068]'
ncfile.notes= "CARM --> Flow measured in a time point that date \n" \
"CHS --> Flow measured each 5' by CHS authomatic gauge; daily mean of flows \n" \
"The discharges are calculated integrating the instantaneous flow for the whole day"
# Crear dimensiones
ncfile.createDimension('date', len(data))
ncfile.createDimension('dif_disch', len(data))
for dim in ncfile.dimensions.items():
    print(dim)
# %%
date_var = ncfile.createVariable('date', np.float32, ('date'))
date_var.long_name = 'Tiempo en días desde 1 enero 1970'
date_var.units = 'day'
date_var[:] = dias_desde_1970.values

flow_var = ncfile.createVariable('dif_disch', np.float32, ('date',))
flow_var.units = 'm³'
flow_var.long_name = 'Accumulated daily discharge difference CHS and CARM '
flow_var[:] = data["Dif_disch"].values
# %%
ncfile.close()
# %% COMPROBACION
dataset = Dataset('albujon_discharge_daily_dif_disch.nc', "r")
print(dataset)
print(dataset.variables.keys())  # Ver las variables en el archivo
print('los fechas data son')
print(dataset.variables["date"][:])
print('los flow  son ')
print(dataset.variables['dif_disch'][:])

# %% Leer las variables
tiempo = dataset.variables["date"][:]  # Días desde 1970
flujo = dataset.variables["dif_disch"][:]    # Flujo en m3
# Convertir tiempo a formato datetime para mejor visualización
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()
# Hacer el plot
plt.figure(figsize=(10, 5))
plt.plot(fechas, flujo, marker="o", linestyle="-", color="b", label="Dif Discharge (m³)")
plt.xlabel("Fecha");plt.ylabel("Dif Discharge (m³)");plt.title("Serie Temporal del dif discharge ");plt.xticks(rotation=45);plt.legend();plt.grid();plt.show()
# %%
