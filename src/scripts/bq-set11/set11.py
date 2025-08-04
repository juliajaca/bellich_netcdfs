# %%
from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import pandas as pd
import shutil
import matplotlib.pyplot as plt
import sys  
sys.path.append("../")
from generar_txt import generar_txt

# %%
data = pd.read_excel('../../datos/bq/set11/Datos_historicos_BELICH.xlsx', )
# %%
# Mapeo de nombres de meses a números
meses_dict = {
    "enero": 1, "febrero": 2, "marzo": 3, "abril": 4, "mayo": 5, "jnio": 6, "junio": 6,
    "julio": 7, "agosto": 8, "sitiembre": 9, "setiembre":9, "octubre": 10, "noviembre": 11, "deciembre": 12, 'diciembre':12,
    "ENERO": 1, "FEBRERO": 2, "MARZO": 3, "ABRIL": 4, "MAYO": 5, "JUNIO": 6,
    "JULIO": 7, "AGOSTO": 8, "SETIEMBRE": 9, "OCTUBRE": 10, "NOVIEMBRE": 11, "DICIEMBRE": 12
}

# Convertir a fecha en el día 15 de cada mes
data["Fecha"] = pd.to_datetime(data["Unnamed: 0"].astype(str) + "-" + data["Unnamed: 1"].str.strip().map(meses_dict).astype(str) + "-1")

data = data.sort_values(by=["Fecha"]).reset_index(drop=True)
data["Fecha"] = pd.to_datetime(data["Fecha"])
data["Fecha_fin"] = data["Fecha"] + pd.offsets.MonthEnd(0)
data = data.replace(np.nan, -9999) #reemplazo los nan por -9999 

# Fecha de referencia (1 de enero de 1970)
epoch = pd.Timestamp("1970-01-01")

data["Días desde 1970"] = (data["Fecha"] - epoch) / pd.Timedelta(days=1)

start_days = (data["Fecha"] - epoch) / pd.Timedelta(days=1)
end_days = (data["Fecha_fin"] - epoch) / pd.Timedelta(days=1)

# Mostrar el resultado
print(data)
# %%
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
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/chlorophyll/BELICH_BIOGQ_V3/'
nombre_fichero = 'BELICH_BIOGQ_V3_CHL'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

ncfile.title= nombre_fichero
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'BELICH_BIOGQ_v3'
ncfile.project = 'BELICH'
ncfile.source = 'In situ data collection'
ncfile.Conventions = 'CF-1.8'
ncfile.comments= "Monthly average for the entire lagoon"

# Crear dimensiones
ncfile.createDimension('time', len(data))
ncfile.createDimension('nv', 2)  # para los bounds

for dim in ncfile.dimensions.items():
    print(dim)
# %%
date_var = ncfile.createVariable('time', np.float64, ('time'))
date_var.units= "days since 1970-01-01 00:00:0"
date_var.calendar = 'gregorian'
date_var.standard_name = "time"
date_var.bounds = "time_bnds"
date_var[:] = data['Días desde 1970'].values

chl_var = ncfile.createVariable('chlorophyll', np.float64, ('time',))
chl_var.units = 'ug L-1'
chl_var.standard_name = 'mass_concentration_of_chlorophyll_a_in_sea_water'
chl_var.long_name = 'Chlorophyll-a Concentration in Sea Water'
chl_var.cell_methods = 'time: mean'
chl_var.missing_value = -9999
chl_var[:] = data["Clorofila (microg L-1) Promedio mensual para toda la laguna"].values

# chl_sd_var = ncfile.createVariable('chl_sd', np.float64, ('time',))
# chl_sd_var.units = 'ug L-1'
# chl_sd_var.standard_name = 'mass_concentration_of_chlorophyll_a_in_sea_water'
# chl_sd_var.cell_methods = 'time: standard_deviation'
# chl_sd_var[:] = data["Clorofila  (microg L-1)  SD"].values

time_bnds = ncfile.createVariable("time_bnds", np.float64, ("time", "nv"))
time_bnds[:, 0] = start_days.values
time_bnds[:, 1] = end_days.values
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
chl = dataset.variables["chlorophyll"][:] 
# chl_sd = dataset.variables["chl_sd"][:] 
bounds = dataset.variables['time_bnds'][:,:] 

print(tiempo)
print('-----------------')
print(chl)
print('-----------------')
# print(chl_sd)
print('-----------------')
print(bounds)

# Convertir tiempo a formato datetime para mejor visualización
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()

# %%
# Hacer los plots
# orden = np.argsort(fechas)  # Obtiene los índices que ordenarían 'fechas'
# fechas_ordenadas = np.array(fechas)[orden]
# chl_ordenado = np.array(chl)[orden]
# chl_sd_ordenado = np.array(chl_sd)[orden]
plt.figure(figsize=(10, 5))
plt.plot(fechas, chl, marker="o", linestyle="-", color="b", label="mean (microg/L)")
plt.xlabel("Fecha"); plt.ylabel("mean chl (microg/L)"); plt.title("Serie Temporal de la clorofila"); plt.xticks(rotation=45); plt.legend(); plt.grid(); plt.show()
# %%
# plt.figure(figsize=(10, 5))
# plt.plot(fechas, chl_sd, marker="o", linestyle="-", color="b", label="SD (microg/L)")
# plt.xlabel("Fecha");plt.ylabel("sd chl (microg/L)");plt.title("Serie Temporal del sd de la clorofila");plt.xticks(rotation=45);plt.legend();plt.grid();plt.show()
# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')

# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Biogeochemical/chlorophyll/BELICH_BIOGQ_V3/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')
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
#  NITRATO
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/nutrients/nitrate/BELICH_BIOGQ_V3/'

nombre_fichero = 'BELICH_BIOGQ_V3_NO3'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

ncfile.title= nombre_fichero
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'BELICH_BIOGQ_v3'
ncfile.project = 'BELICH'
ncfile.source = 'In situ data collection'
ncfile.Conventions = 'CF-1.8'
ncfile.comments= "Monthly average for the entire lagoon"

# Crear dimensiones
ncfile.createDimension('time', len(data))
ncfile.createDimension('nv', 2)  # para los bounds

for dim in ncfile.dimensions.items():
    print(dim)
# %%
date_var = ncfile.createVariable('time', np.float64, ('time'))
date_var.units= "days since 1970-01-01 00:00:0"
date_var.calendar = 'gregorian'
date_var.standard_name = "time"
date_var.bounds = "time_bnds"
date_var[:] = data['Días desde 1970'].values

nitrato_var = ncfile.createVariable('nitrate', np.float64, ('time',))
nitrato_var.units = 'umol L-1'
nitrato_var.standard_name = 'mole_concentration_of_nitrate_in_sea_water'
nitrato_var.long_name = 'Nitrate concentration in sea water'
nitrato_var.cell_methods = 'time: mean'
nitrato_var.missing_value = -9999
nitrato_var[:] = data["Nitrato (microM) Promedio"].values

# nitrato_sd_var = ncfile.createVariable('nitrate_sd', np.float64, ('time',))
# nitrato_sd_var.units = 'umol L-1'
# nitrato_sd_var.cell_methods = 'time: standard_deviation'
# nitrato_sd_var.standard_name = 'mole_concentration_of_nitrate_in_sea_water'
# nitrato_sd_var[:] = data["Nitrato (microM) SD"].values

time_bnds = ncfile.createVariable("time_bnds", np.float64, ("time", "nv"))
time_bnds[:, 0] = start_days.values
time_bnds[:, 1] = end_days.values
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
chl = dataset.variables["nitrate"][:] 
# chl_sd = dataset.variables["nitrate_sd"][:]    # 

print(tiempo)
print('-----------------')
print(chl)
print('-----------------')
# print(chl_sd)

# Convertir tiempo a formato datetime para mejor visualización
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()

# %%
# Hacer los plots
orden = np.argsort(fechas)  # Obtiene los índices que ordenarían 'fechas'
fechas_ordenadas = np.array(fechas)[orden]
chl_ordenado = np.array(chl)[orden]
# chl_sd_ordenado = np.array(chl_sd)[orden]
plt.figure(figsize=(10, 5))
plt.plot(fechas_ordenadas, chl_ordenado, marker="o", linestyle="-", color="b", label="mean (microM/L)")
plt.xlabel("Fecha"); plt.ylabel("mean nitrato (microM/L)"); plt.title("Serie Temporal del nitrato"); plt.xticks(rotation=45); plt.legend(); plt.grid(); plt.show()
# %%
# plt.figure(figsize=(10, 5))
# plt.plot(fechas_ordenadas, chl_sd_ordenado, marker="o", linestyle="-", color="b", label="SD (microM/L)")
# plt.xlabel("Fecha");plt.ylabel("sd nitrato (microM/L)");plt.title("Serie Temporal del sd del nitrato");plt.xticks(rotation=45);plt.legend();plt.grid();plt.show()
# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')
# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Biogeochemical/nutrients/nitrate/BELICH_BIOGQ_V3/'

shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')
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
# %% NITRITO
nombre_fichero = 'BELICH_BIOGQ_V3_NO2'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/nutrients/nitrite/BELICH_BIOGQ_V3/'

ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

ncfile.title= nombre_fichero
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'BELICH_BIOGQ_v3'
ncfile.project = 'BELICH'
ncfile.source = 'In situ data collection'
ncfile.Conventions = 'CF-1.8'
ncfile.comments= "Monthly average for the entire lagoon"

# Crear dimensiones
ncfile.createDimension('time', len(data))
ncfile.createDimension('nv', 2)  # para los bounds

for dim in ncfile.dimensions.items():
    print(dim)
# %%
date_var = ncfile.createVariable('time', np.float64, ('time'))
date_var.units= "days since 1970-01-01 00:00:0"
date_var.calendar = 'gregorian'
date_var.standard_name = "time"
date_var.bounds = "time_bnds"
date_var[:] = data['Días desde 1970'].values

nitrato_var = ncfile.createVariable('nitrite', np.float64, ('time',))
nitrato_var.units = 'umol L-1'
nitrato_var.standard_name = 'mole_concentration_of_nitrite_in_sea_water'
nitrato_var.long_name = 'Nitrite concentration in sea water'
nitrato_var.cell_methods = 'time: mean'
nitrato_var.missing_value = -9999
nitrato_var[:] = data["Nitrito (microM) Promedio"].values

# nitrato_sd_var = ncfile.createVariable('nitrite_sd', np.float64, ('time',))
# nitrato_sd_var.units = 'umol L-1'
# nitrato_sd_var.cell_methods = 'time: standard_deviation'
# nitrato_sd_var.standard_name = 'mole_concentration_of_nitrite_in_sea_water'
# nitrato_sd_var[:] = data["Nitrito (microM) SD"].values

time_bnds = ncfile.createVariable("time_bnds", np.float64, ("time", "nv"))
time_bnds[:, 0] = start_days.values
time_bnds[:, 1] = end_days.values
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
chl = dataset.variables["nitrite"][:] 
# chl_sd = dataset.variables["nitrite_sd"][:]    # 

print(tiempo)
print('-----------------')
print(chl)
print('-----------------')
# print(chl_sd)

# Convertir tiempo a formato datetime para mejor visualización
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()

# %%
# Hacer los plots
orden = np.argsort(fechas)  # Obtiene los índices que ordenarían 'fechas'
fechas_ordenadas = np.array(fechas)[orden]
chl_ordenado = np.array(chl)[orden]
# chl_sd_ordenado = np.array(chl_sd)[orden]
plt.figure(figsize=(10, 5))
plt.plot(fechas_ordenadas, chl_ordenado, marker="o", linestyle="-", color="b", label="mean (microM/L)")
plt.xlabel("Fecha"); plt.ylabel("mean nitrito (microM/L)"); plt.title("Serie Temporal del nitrito"); plt.xticks(rotation=45); plt.legend(); plt.grid(); plt.show()
# %%
# plt.figure(figsize=(10, 5))
# plt.plot(fechas_ordenadas, chl_sd_ordenado, marker="o", linestyle="-", color="b", label="SD (microM/L)")
# plt.xlabel("Fecha");plt.ylabel("sd nitrito (microM/L)");plt.title("Serie Temporal del sd del nitrito");plt.xticks(rotation=45);plt.legend();plt.grid();plt.show()
# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')
# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Biogeochemical/nutrients/nitrite/BELICH_BIOGQ_V3/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')
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
"""
# %% FOSFATO
nombre_fichero = 'BELICH_BIOGQ_V3_PO4'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/nutrients/phosphate/BELICH_BIOGQ_V3/'

ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

ncfile.title= nombre_fichero
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'BELICH_BIOGQ_v3'
ncfile.project = 'BELICH'
ncfile.source = 'In situ data collection'
ncfile.Conventions = 'CF-1.8'
ncfile.comments= "Monthly average for the entire lagoon"

# Crear dimensiones
ncfile.createDimension('time', len(data))
ncfile.createDimension('nv', 2)  # para los bounds
for dim in ncfile.dimensions.items():
    print(dim)
# %%
date_var = ncfile.createVariable('time', np.float64, ('time'))
date_var.units= "days since 1970-01-01 00:00:0"
date_var.calendar = 'gregorian'
date_var.standard_name = "time"
date_var.bounds = "time_bnds"
date_var[:] = data['Días desde 1970'].values

nitrato_var = ncfile.createVariable('phosphate', np.float64, ('time',))
nitrato_var.units = 'umol L-1'
nitrato_var.standard_name = 'mole_concentration_of_phosphate_in_sea_water'
nitrato_var.long_name = 'Phosphate concentration in sea water'
nitrato_var.cell_methods = 'time: mean'
nitrato_var.missing_value = -9999
nitrato_var[:] = data["Fosfato (microM) Promedio"].values

# nitrato_sd_var = ncfile.createVariable('phosphate_sd', np.float64, ('time',))
# nitrato_sd_var.units = 'umol L-1'
# nitrato_sd_var.cell_methods = 'time: standard_deviation'
# nitrato_sd_var.standard_name = 'mole_concentration_of_phosphate_in_sea_water'
# nitrato_sd_var[:] = data["Fosfato (microM) SD"].values

time_bnds = ncfile.createVariable("time_bnds", np.float64, ("time", "nv"))
time_bnds[:, 0] = start_days.values
time_bnds[:, 1] = end_days.values
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
chl = dataset.variables["phosphate"][:] 
# chl_sd = dataset.variables["phosphate_sd"][:]    # 

print(tiempo)
print('-----------------')
print(chl)
print('-----------------')
# print(chl_sd)

# Convertir tiempo a formato datetime para mejor visualización
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()

# %%
# Hacer los plots
orden = np.argsort(fechas)  # Obtiene los índices que ordenarían 'fechas'
fechas_ordenadas = np.array(fechas)[orden]
chl_ordenado = np.array(chl)[orden]
# chl_sd_ordenado = np.array(chl_sd)[orden]
plt.figure(figsize=(10, 5))
plt.plot(fechas_ordenadas, chl_ordenado, marker="o", linestyle="-", color="b", label="mean (microM/L)")
plt.xlabel("Fecha"); plt.ylabel("mean fosfato (microM/L)"); plt.title("Serie Temporal del fosfato"); plt.xticks(rotation=45); plt.legend(); plt.grid(); plt.show()
# %%
# plt.figure(figsize=(10, 5))
# plt.plot(fechas_ordenadas, chl_sd_ordenado, marker="o", linestyle="-", color="b", label="SD (microM/L)")
# plt.xlabel("Fecha");plt.ylabel("sd fosfato (microM/L)");plt.title("Serie Temporal del sd del fosfato");plt.xticks(rotation=45);plt.legend();plt.grid();plt.show()
# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')
# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Biogeochemical/nutrients/phosphate/BELICH_BIOGQ_V3/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')
# %%
