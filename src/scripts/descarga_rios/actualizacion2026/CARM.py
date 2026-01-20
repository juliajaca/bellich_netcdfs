# %%
from netCDF4 import Dataset, stringtochar
import numpy as np
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("../../")
from generar_txt import generar_txt
# %%
n_dataset ='CARM_Albujon'
nombre_fichero = 'CARM_PR'
data = pd.read_csv('C:/Users/Julia/Documents/VSCODE_BELLICH/src/scripts/descarga_rios/actualizacion2026/datos/CARM/CARM_Albujon_flow_water_chemistry.csv',  sep=';', parse_dates= ['Date'],) #son todo strings
# data.groupby(['Variable']).size()

"""
 ******** **         *******   **       **
/**///// /**        **/////** /**      /**
/**      /**       **     //**/**   *  /**
/******* /**      /**      /**/**  *** /**
/**////  /**      /**      /**/** **/**/**
/**      /**      //**     ** /**** //****
/**      /******** //*******  /**/   ///**
//       ////////   ///////   //       //
"""
data_flow = data[['Date', 'Flow']]
data_flow = data_flow.dropna(subset=['Flow'])

# %%
nombre_fichero = 'CARM_Albujon_FLOW'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Runoff/CARM_Albujon/'
path_copia = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Runoff/'

ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')

ncfile.title = nombre_fichero
ncfile.institution="CARM"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = n_dataset
ncfile.project = n_dataset
ncfile.source = 'In situ data collection'
ncfile.Conventions = 'CF-1.8'
ncfile.summary = 'Time series of water flow and water chemistry variables measured at the Rambla del Albujón outlet (Mar Menor watershed, SE Spain). Data include instantaneous flow, electrical conductivity, nitrate and phosphate concentrations, collected between August 2019 and December 2025.'
ncfile.comment ='Sampling was carried out in the last 200 m of the Rambla del Albujón watercourse, between latitudes 37.716221-37.716013 and longitudes -0.861232--0.858857 (WGS84).\n Sampling interval was approximately weekly until August 2021 and every 2-3 days thereafter. Sampling hour was not recorded but usually corresponds to the first light hours of the day.\n In some sampling campaigns multiple localities were sampled over two or three consecutive days. In these cases, the sampling date assigned to the Rambla del Albujón outlet corresponds to the first day of the campaign, as this location is usually sampled first. \nMissing values are coded as NA in the original data source.'

# Variable time
data_flow = data_flow.sort_values(by=["Date"]).reset_index(drop=True)

unique_times = data_flow['Date'].drop_duplicates().sort_values().reset_index(drop=True)
fechas_unicas_ts = pd.to_datetime(unique_times)

# ANALISIS DUPLICADOS
# Número total de filas
print("Filas totales:", len(data_flow))
# Número de fechas únicas
print("Fechas únicas:", data_flow['Date'].nunique())
duplicadas = data_flow[data_flow.duplicated(subset='Date', keep=False)]
print(duplicadas.sort_values('Date'))
# -----------------------------------------------------------

epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (fechas_unicas_ts - epoch) / pd.Timedelta(days=1)

ncfile.createDimension('time', len(dias_desde_1970))
for dim in ncfile.dimensions.items():
    print(dim)
# %%
time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:00"
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values # Se asigna directamente

valores_con_nan = data_flow.copy()
valores_con_nan[np.isnan(valores_con_nan)] = -9999
value_var = ncfile.createVariable('flow', np.float64, ('time', ))
value_var.units= 'L s-1'
value_var.standard_name = "water_volume_transport_into_sea_water"
value_var.long_name = 'Instantaneous water flow'
value_var.missing_value = -9999
value_var[:] = valores_con_nan['Flow']
value_var.comment = 'Flow was manually measured by delineating the channel cross-section and measuring current velocity using a current meter (instrument model unknown)'

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

tiempo = dataset.variables["time"][:]  # Días desde 1970
value = dataset.variables['flow'][:] 

print(tiempo); print('-----------------')
print(value[0:10]); print('-----------------')

fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()

# %% PLOT
fechas = np.array(fechas)  # convertir a numpy array para indexar
orden = np.argsort(fechas)  # obtener índices ordenados

# Reordenar fechas 
fechas_ordenadas = fechas[orden]

fig, axes = plt.subplots(1, 1, figsize=(40,20), sharex=True, sharey=True)
axes.plot(fechas, value, marker='o', linestyle='None',markersize=4 )
axes.grid(True)
axes.set_ylabel("instant flow")
axes.set_xlabel('Fecha')
plt.tight_layout()
plt.savefig(
    f"{path}CARM_flow.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()

# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')
# %%
shutil.copy(f'{path}{nombre_fichero}.nc',f'{path_copia}{nombre_fichero}.nc')
# %%
"""
. . .-. .-. .-. .-. .-. .-.
|\|  |   |  |(  |-|  |  |-
' ` `-'  '  ' ' ` '  '  `-'

"""
data_flow = data[['Date', 'Nitrate']]
data_flow = data_flow.dropna(subset=['Nitrate'])

# %%
nombre_fichero = 'CARM_Albujon_NO3'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Runoff/CARM_Albujon/'
path_copia = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Runoff/'

ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')

ncfile.title = nombre_fichero
ncfile.institution = "CARM"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = n_dataset
ncfile.project = n_dataset
ncfile.source = 'In situ data collection'
ncfile.Conventions = 'CF-1.8'
ncfile.summary = 'Time series of water flow and water chemistry variables measured at the Rambla del Albujón outlet (Mar Menor watershed, SE Spain). Data include instantaneous flow, electrical conductivity, nitrate and phosphate concentrations, collected between August 2019 and December 2025.'
ncfile.comment ='Sampling was carried out in the last 200 m of the Rambla del Albujón watercourse, between latitudes 37.71622-37.716013 and longitudes -0.861232--0.858857 (WGS84).\n Sampling interval was approximately weekly until August 2021 and every 2-3 days thereafter. Sampling hour was not recorded but usually corresponds to the first light hours of the day.\n In some sampling campaigns multiple localities were sampled over two or three consecutive days. In these cases, the sampling date assigned to the Rambla del Albujón outlet corresponds to the first day of the campaign, as this location is usually sampled first. \nMissing values are coded as NA in the original data source.'

# Variable time
data_flow = data_flow.sort_values(by=["Date"]).reset_index(drop=True)

unique_times = data_flow['Date'].drop_duplicates().sort_values().reset_index(drop=True)
fechas_unicas_ts = pd.to_datetime(unique_times)

# ANALISIS DUPLICADOS
# Número total de filas
print("Filas totales:", len(data_flow))
# Número de fechas únicas
print("Fechas únicas:", data_flow['Date'].nunique())
duplicadas = data_flow[data_flow.duplicated(subset='Date', keep=False)]
print(duplicadas.sort_values('Date'))
# -----------------------------------------------------------

epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (fechas_unicas_ts - epoch) / pd.Timedelta(days=1)

ncfile.createDimension('time', len(dias_desde_1970))
for dim in ncfile.dimensions.items():
    print(dim)
# %%
time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:00"
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values # Se asigna directamente

valores_con_nan = data_flow.copy()
valores_con_nan[np.isnan(valores_con_nan)] = -9999
value_var = ncfile.createVariable('nitrate', np.float64, ('time', ))
value_var.units= 'mg L-1'
value_var.standard_name = 'mass_concentration_of_nitrate_in_sea_water'
value_var.long_name = 'Dissolved nitrate concentration in sea water'
value_var.missing_value = -9999
value_var[:] = valores_con_nan['Nitrate']
value_var.comment = 'Nitrate concentration was measured in the laboratory. Analytical methods are not reported in the original metadata.'

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

tiempo = dataset.variables["time"][:]  # Días desde 1970
value = dataset.variables['nitrate'][:] 

print(tiempo); print('-----------------')
print(value[0:10]); print('-----------------')

fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()

# %% PLOT
fechas = np.array(fechas)  # convertir a numpy array para indexar
orden = np.argsort(fechas)  # obtener índices ordenados

# Reordenar fechas 
fechas_ordenadas = fechas[orden]

fig, axes = plt.subplots(1, 1, figsize=(40,20), sharex=True, sharey=True)
axes.plot(fechas, value, marker='o', linestyle='None',markersize=4 )
axes.grid(True)
axes.set_ylabel("nitrate")
axes.set_xlabel('Fecha')
plt.tight_layout()
plt.savefig(
    f"{path}CARM_nitrate.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()

# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')
# %%
shutil.copy(f'{path}{nombre_fichero}.nc',f'{path_copia}{nombre_fichero}.nc')


"""
.sSSSSs.    .sSSSSs.    .sSSSSs.    .sSSSSs.    .sSSSSs.       .sSSSSSSSSs.   .sSSSSs.
SSSSSSSSSs. SSSSSSSSSs. SSSSSSSSSs. SSSSSSSSSs. SSSSSSSSSs. .sSSSSSSSSSSSSSs. SSSSSSSSSs.
S SSS SSSS' S SSS SSSSS S SSS SSSS' S SSS SSSS' S SSS SSSSS SSSSS S SSS SSSSS S SSS SSSSS
S  SS       S  SS SSSSS S  SS       S  SS       S  SS SSSSS SSSSS S  SS SSSSS S  SS SSSSS
S..SSsss    S..SS SSSSS `SSSSsSSSa. S..SSsss    S..SSsSSSSS `:S:' S..SS `:S:' S..SS SSSSS
S:::SSSS    S:::S SSSSS .sSSS SSSSS S:::SSSS    S:::S SSSSS       S:::S       S:::S SSSSS
S;;;S       S;;;S SSSSS S;;;S SSSSS S;;;S       S;;;S SSSSS       S;;;S       S;;;S SSSSS
S%%%S       S%%%S SSSSS S%%%S SSSSS S%%%S       S%%%S SSSSS       S%%%S       S%%%S SSSSS
SSSSS       SSSSSsSSSSS SSSSSsSSSSS SSSSS       SSSSS SSSSS       SSSSS       SSSSSsSSSSS

"""
data_flow = data[['Date', 'Phosphate']]
data_flow = data_flow.dropna(subset=['Phosphate'])

# %%
nombre_fichero = 'CARM_Albujon_PO4'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Runoff/CARM_Albujon/'
path_copia = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Runoff/'

ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')

ncfile.title = nombre_fichero
ncfile.institution = "CARM"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = n_dataset
ncfile.project = n_dataset
ncfile.source = 'In situ data collection'
ncfile.Conventions = 'CF-1.8'
ncfile.summary = 'Time series of water flow and water chemistry variables measured at the Rambla del Albujón outlet (Mar Menor watershed, SE Spain). Data include instantaneous flow, electrical conductivity, nitrate and phosphate concentrations, collected between August 2019 and December 2025.'
ncfile.comment ='Sampling was carried out in the last 200 m of the Rambla del Albujón watercourse, between latitudes 37.716221-37.716013 and longitudes -0.861232--0.858857 (WGS84).\n Sampling interval was approximately weekly until August 2021 and every 2-3 days thereafter. Sampling hour was not recorded but usually corresponds to the first light hours of the day.\n In some sampling campaigns multiple localities were sampled over two or three consecutive days. In these cases, the sampling date assigned to the Rambla del Albujón outlet corresponds to the first day of the campaign, as this location is usually sampled first. \nMissing values are coded as NA in the original data source.'

# Variable time
data_flow = data_flow.sort_values(by=["Date"]).reset_index(drop=True)

unique_times = data_flow['Date'].drop_duplicates().sort_values().reset_index(drop=True)
fechas_unicas_ts = pd.to_datetime(unique_times)

# ANALISIS DUPLICADOS
# Número total de filas
print("Filas totales:", len(data_flow))
# Número de fechas únicas
print("Fechas únicas:", data_flow['Date'].nunique())
duplicadas = data_flow[data_flow.duplicated(subset='Date', keep=False)]
print(duplicadas.sort_values('Date'))
# -----------------------------------------------------------

epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (fechas_unicas_ts - epoch) / pd.Timedelta(days=1)

ncfile.createDimension('time', len(dias_desde_1970))
for dim in ncfile.dimensions.items():
    print(dim)
# %%
time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:00"
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values # Se asigna directamente

valores_con_nan = data_flow.copy()
valores_con_nan[np.isnan(valores_con_nan)] = -9999
value_var = ncfile.createVariable('phosphate', np.float64, ('time', ))
value_var.units= 'mg L-1'
value_var.standard_name = 'mass_concentration_of_phosphate_in_sea_water'
value_var.long_name = 'Dissolved phosphate concentration in sea water'
value_var.missing_value = -9999
value_var[:] = valores_con_nan['Phosphate']
value_var.comment = 'Phosphate concentration was measured in the laboratory (analytical methods not reported). Concentrations below the limit of detection were recoded to zero in the original dataset.\nThe reported limit of detection was 0.061 mg PO4 L-1 until 17 April 2024, and 0.2 mg PO4 L-1 from that date onwards.\nPhosphate data are available from 23 August 2021 onwards.'

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

tiempo = dataset.variables["time"][:]  # Días desde 1970
value = dataset.variables['phosphate'][:] 

print(tiempo); print('-----------------')
print(value[0:10]); print('-----------------')

fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()

# %% PLOT
fechas = np.array(fechas)  # convertir a numpy array para indexar
orden = np.argsort(fechas)  # obtener índices ordenados

# Reordenar fechas 
fechas_ordenadas = fechas[orden]

fig, axes = plt.subplots(1, 1, figsize=(40,20), sharex=True, sharey=True)
axes.plot(fechas, value, marker='o', linestyle='None',markersize=4 )
axes.grid(True)
axes.set_ylabel("fosfato")
axes.set_xlabel('Fecha')
plt.tight_layout()
plt.savefig(
    f"{path}CARM_fosfato.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()

# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')
# %%
shutil.copy(f'{path}{nombre_fichero}.nc',f'{path_copia}{nombre_fichero}.nc')

"""
.sSSSSs.    .sSSSSs.    .sSSSs.  SSSSS .sSSSSs.    .sSSS s.    .sSSSSs.       .sSSSSSSSSs.
SSSSSSSSSs. SSSSSSSSSs. SSSSS SS SSSSS SSSSSSSSSs. SSSSS SSSs. SSSSSSSSSs. .sSSSSSSSSSSSSSs.
S SSS SSSSS S SSS SSSSS S SSS  `sSSSSS S SSS SSSSS S SSS SSSSS S SSS SSSSS SSSSS S SSS SSSSS
S  SS SSSS' S  SS SSSSS S  SS    SSSSS S  SS SSSSS S  SS SSSSS S  SS SSSS' SSSSS S  SS SSSSS
S..SS       S..SS SSSSS S..SS    SSSSS S..SS SSSSS S..SS SSSSS S..SS       `:S:' S..SS `:S:'
S:::S SSSSS S:::S SSSSS S:::S    SSSSS S:::S SSSSS S:::S SSSSS S:::S SSSSS       S:::S
S;;;S SSSSS S;;;S SSSSS S;;;S    SSSSS S;;;S SSSSS S;;;S SSSSS S;;;S SSSSS       S;;;S
S%%%S SSSSS S%%%S SSSSS S%%%S    SSSSS S%%%S SSSS' S%%%S SSSSS S%%%S SSSSS       S%%%S
SSSSSsSSSSS SSSSSsSSSSS SSSSS    SSSSS SSSSSsS;:'  SSSSSsSSSSS SSSSSsSSSSS       SSSSS

"""
data_flow = data[['Date', 'EConductivity']]
data_flow = data_flow.dropna(subset=['EConductivity'])

# %%
nombre_fichero = 'CARM_Albujon_COND'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Runoff/CARM_Albujon/'
path_copia = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Runoff/'

ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')

ncfile.title = nombre_fichero
ncfile.institution = "CARM"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = n_dataset
ncfile.project = n_dataset
ncfile.source = 'In situ data collection'
ncfile.Conventions = 'CF-1.8'
ncfile.summary = 'Time series of water flow and water chemistry variables measured at the Rambla del Albujón outlet (Mar Menor watershed, SE Spain). Data include instantaneous flow, electrical conductivity, nitrate and phosphate concentrations, collected between August 2019 and December 2025.'
ncfile.comment ='Sampling was carried out in the last 200 m of the Rambla del Albujón watercourse, between latitudes 37.716221-37.716013 and longitudes -0.861232--0.858857 (WGS84).\n Sampling interval was approximately weekly until August 2021 and every 2-3 days thereafter. Sampling hour was not recorded but usually corresponds to the first light hours of the day.\n In some sampling campaigns multiple localities were sampled over two or three consecutive days. In these cases, the sampling date assigned to the Rambla del Albujón outlet corresponds to the first day of the campaign, as this location is usually sampled first. \nMissing values are coded as NA in the original data source.'

# Variable time
data_flow = data_flow.sort_values(by=["Date"]).reset_index(drop=True)

unique_times = data_flow['Date'].drop_duplicates().sort_values().reset_index(drop=True)
fechas_unicas_ts = pd.to_datetime(unique_times)

# ANALISIS DUPLICADOS
# Número total de filas
print("Filas totales:", len(data_flow))
# Número de fechas únicas
print("Fechas únicas:", data_flow['Date'].nunique())
duplicadas = data_flow[data_flow.duplicated(subset='Date', keep=False)]
print(duplicadas.sort_values('Date'))
# -----------------------------------------------------------

epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (fechas_unicas_ts - epoch) / pd.Timedelta(days=1)

ncfile.createDimension('time', len(dias_desde_1970))
for dim in ncfile.dimensions.items():
    print(dim)
# %%
time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:00"
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values # Se asigna directamente

valores_con_nan = data_flow.copy()
valores_con_nan[np.isnan(valores_con_nan)] = -9999
value_var = ncfile.createVariable('conductivity', np.float64, ('time', ))
value_var.units= 'µS cm-1'
value_var.standard_name = "sea_water_electrical_conductivity"
value_var.long_name = "Electrical conductivity of water"
value_var.missing_value = -9999
value_var[:] = valores_con_nan['EConductivity']
value_var.comment = 'Electrical conductivity was measured in situ using a conductivity meter (instrument model unknown). Values are assumed to be temperature-corrected to 20 °C, although this could not be confirmed.'

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

tiempo = dataset.variables["time"][:]  # Días desde 1970
value = dataset.variables['conductivity'][:] 

print(tiempo); print('-----------------')
print(value[0:10]); print('-----------------')

fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()

# %% PLOT
fechas = np.array(fechas)  # convertir a numpy array para indexar
orden = np.argsort(fechas)  # obtener índices ordenados

# Reordenar fechas 
fechas_ordenadas = fechas[orden]

fig, axes = plt.subplots(1, 1, figsize=(40,20), sharex=True, sharey=True)
axes.plot(fechas, value, marker='o', linestyle='None',markersize=4 )
axes.grid(True)
axes.set_ylabel("CONDUCTIVIDAD")
axes.set_xlabel('Fecha')
plt.tight_layout()
plt.savefig(
    f"{path}CARM_COND.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()
# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}_display.txt')
# %%
shutil.copy(f'{path}{nombre_fichero}.nc',f'{path_copia}{nombre_fichero}.nc')
# %%