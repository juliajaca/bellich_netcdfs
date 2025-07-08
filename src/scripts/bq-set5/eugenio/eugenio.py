# %%
from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from geopy.distance import geodesic
import sys
sys.path.append("../../")
from generar_txt import generar_txt
# %%
"""
     **     *******     ********   ****   **   ****   ****
    ****   /**////**   **//////** */// * ***  */// * */// *
   **//**  /**   /**  **      // /*   /*//** /*   /*/    /*
  **  //** /*******  /**         / ****  /** / ****    ***
 **********/**///**  /**    ***** */// * /**  */// *  *//
/**//////**/**  //** //**  ////**/*   /* /** /*   /* *
/**     /**/**   //** //******** / ****  ****/ **** /******
//      // //     //   ////////   ////  ////  ////  //////
"""
'''
08/1981  he puesto 01/08/1981
'''
data0 = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/bq/set5/MEDIAS clorofila_histórico para Eugenio.xlsx',  dtype={  
"A": "float64", "B": "float64", "C": "float64", }, parse_dates= ['FECHA'], sheet_name = 'ARG8182', skiprows=1,  nrows= 24, usecols= ['FECHA', 'modo', 'A', 'B', 'C', ])
print(data0.head())
print(f'tiene una longitud de {len(data0)} filas')
print(data0.tail())
# %%
data = data0.sort_values("FECHA").reset_index(drop=True)
# %%
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (data.FECHA.drop_duplicates() - epoch) / pd.Timedelta(days=1)
# %%
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/set5/chl eugenio/'

ncfile = Dataset(f'{path}clorofila_eugenio_ARG8182.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)
# 
ncfile.title='clorofila eugenio pagina ARG8182'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'; ncfile.source = 'XXX'; ncfile.Conventions = 'CF-1.8'

# %%
max_length_param = len("A")
estaciones = ['A', 'B', 'C',]
estaciones_np = np.array(estaciones, dtype=f'S{max(len(s) for s in estaciones)}')


df_chl =  data.iloc[:, [2,3,4, ]]
ncfile.createDimension('time', len(dias_desde_1970))
ncfile.createDimension('unit_char_len', max_length_param)
ncfile.createDimension('station_name', len(estaciones))
ncfile.createDimension('depth', 2)

for dim in ncfile.dimensions.items():
    print(dim)

# %%
time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
time_var[:] = dias_desde_1970.values  # Se asigna directamente

parameter_var = ncfile.createVariable('station', 'S1', ('station_name', 'unit_char_len'))
parameter_var.long_name = 'station'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)

mode_var = ncfile.createVariable('depth', np.int8, ('depth'))
mode_var.long_name = 'depth'
mode_var.units = 'm'
mode_var.comment = '0m is surface depth, and -9m is bottom depth'
mode_var[:] =[0,-9]
# %%
value_var = ncfile.createVariable('clorofila', np.float64, ('time','depth', 'station_name'))
value_var.long_name = 'mass_concentration_of_chlorophyll_in_sea_water'
value_var.units= '¿mg L-1?'
# %%
multi_index = pd.MultiIndex.from_product([data.FECHA.drop_duplicates(), [0,-9]], names=["FECHA", "modo"])
pivoted = data.pivot_table(index=["FECHA", "modo"], values=estaciones)
pivoted = pivoted.reindex(multi_index)
chlorophyll_data = pivoted.to_numpy().reshape(12, 2, 3)
# %%
value_var[:,:, :] = chlorophyll_data

lat_var = ncfile.createVariable('latitude', np.float64,('station_name') )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var[:] = [ 9999, 9999, 9999]

lon_var = ncfile.createVariable('longitude', np.float64,('station_name') )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var[:] = [ 9999, 9999, 9999]

crs = ncfile.createVariable("crs", "i") 
crs.long_name=  "Reference coordinate system"
crs.grid_mapping_name = "latitude_longitude"
crs.epsg_code = "EPSG:4326"

# %%
ncfile.close()
# %% COMPROBACION
dataset = Dataset(f'{path}clorofila_eugenio_ARG8182.nc', "r")
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
unit = dataset.variables["station"][:]    #  
value = dataset.variables["clorofila"][:] 
modo = dataset.variables['depth'][:]
print(tiempo); print('-----------------')
print(unit); print('-----------------')
print(value); print('-----------------')
print(modo); print('-----------------')

# %%
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()
# %% PLOT
fechas = np.array(fechas)  # convertir a numpy array para indexar
orden = np.argsort(fechas)  # obtener índices ordenados

# Reordenar fechas y nitrato
fechas_ordenadas = fechas[orden]
chl_ordenado = value[orden, :,:]

fig, axs = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

for i, ax in enumerate(axs):
    ax.plot(fechas_ordenadas, chl_ordenado[:, 0, i], label='S', color='royalblue', marker='o')
    ax.plot(fechas_ordenadas, chl_ordenado[:, 1, i], label='F', color='seagreen', marker = 's')
    
    ax.set_ylabel("Clorofila ")
    ax.set_title(f"Estación {estaciones[i]}")
    ax.legend()
    ax.grid(True)

axs[-1].set_xlabel("Fecha")
plt.tight_layout()
plt.show()

# %%
generar_txt(f'{path}clorofila_eugenio_ARG8182.nc', f'{path}clorofila_eugenio_ARG8182_display.txt')
"""
      **   ******   ****   ****   ****   ****
     /**  **////** */// * *///** */// * */// *
     /** **    // /*   /*/*  */*/*   /*/    /*
     /**/**       / **** /* * /*/ ****    ***
     /**/**        ///*  /**  /* ///*    *//
 **  /**//**    **   *   /*   /*   *    *
//*****  //******   *    / ****   *    /******
 /////    //////   /      ////   /     //////
"""
# %%
'''
Hay un dia que pone 1900 y otro de 1982. Supongo que es un error, cual es la fecha que va en su lugar?
1900 = 1990 y 1982 = 1992
'''
data0 = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/bq/set5/MEDIAS clorofila_histórico para Eugenio.xlsx',  dtype={  
"E1": "float64", "E2": "float64", "E3": "float64","E4": "float64", }, parse_dates= ['fecha'], sheet_name = 'JC9092', skiprows=1,  nrows= 120, usecols= ['fecha', 'E1', 'E2', 'E3', 'E4']) 
print(data0.head())
print(f'tiene una longitud de {len(data0)} filas')
print(data0.tail())
# %%
data = data0.sort_values("fecha").reset_index(drop=True)
# %%
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (data.fecha - epoch) / pd.Timedelta(days=1)
# %%
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/set5/chl eugenio/'
ncfile = Dataset(f'{path}clorofila_eugenio_JC9092.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='clorofila eugenio pçagina JC9092'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'; ncfile.source = 'XXX'; ncfile.Conventions = 'CF-1.8'

# %%
max_length_param = len("E1")
estaciones = ['E1', 'E2', 'E3', 'E4']
estaciones_np = np.array(estaciones, dtype=f'S{max(len(s) for s in estaciones)}')

df_chl =  data.iloc[:, [1,2,3, 4]]
ncfile.createDimension('time', len(dias_desde_1970))
ncfile.createDimension('unit_char_len', max_length_param)
ncfile.createDimension('station_name', len(estaciones))
for dim in ncfile.dimensions.items():
    print(dim)

# %%
time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
time_var[:] = dias_desde_1970.values  # Se asigna directamente

parameter_var = ncfile.createVariable('station', 'S1', ('station_name', 'unit_char_len'))
parameter_var.long_name = 'station'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)
# %%
value_var = ncfile.createVariable('chlorophyll', np.float64, ('time','station_name'))
value_var.long_name = 'mass_concentration_of_chlorophyll_in_sea_water'
value_var.units= 'mg L-1?¿'
value_var[:,:] = df_chl

lat_var = ncfile.createVariable('latitude', np.float64,('station_name') )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var[:] = [ 9999, 9999, 9999, 9999]

lon_var = ncfile.createVariable('longitude', np.float64,('station_name') )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var[:] = [ 9999, 9999, 9999 , 9999]

crs = ncfile.createVariable("crs", "i") 
crs.long_name=  "Reference coordinate system"
crs.grid_mapping_name = "latitude_longitude"
crs.epsg_code = "EPSG:4326"

# %%
ncfile.close()
# %% COMPROBACION
dataset = Dataset(f'{path}clorofila_eugenio_JC9092.nc', "r")
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
unit = dataset.variables["station"][:]    #  
value = dataset.variables["chlorophyll"][:] 
print(tiempo); print('-----------------')
print(unit); print('-----------------')
print(value); print('-----------------')
# %%
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()
# %% PLOT
fechas = np.array(fechas)  # convertir a numpy array para indexar
orden = np.argsort(fechas)  # obtener índices ordenados

# Reordenar fechas y nitrato
fechas_ordenadas = fechas[orden]
chl_ordenado = value[orden, :]

fig, axes = plt.subplots(nrows=len(estaciones), figsize=(10, 10), sharex=True)

for i in range(len(estaciones)):
    axes[i].plot(fechas_ordenadas, chl_ordenado[:, i], marker='o', label=estaciones[i])
    axes[i].set_ylabel('clorofila ')
    axes[i].set_title(f'clorofila_eugenio_JC9092- {estaciones[i]}')
    axes[i].grid(True)
    axes[i].legend()

axes[-1].set_xlabel('Fecha')
plt.tight_layout()
plt.show()

# %%
generar_txt(f'{path}clorofila_eugenio_JC9092.nc', f'{path}clorofila_eugenio_JC9092_display.txt')
# %%

"""
 ******** **     ** *******     *******     ********  ******** **
/**///// /**    /**/**////**   **/////**   **//////**/**///// /**
/**      /**    /**/**   /**  **     //** **      // /**      /**
/******* /**    /**/*******  /**      /**/**         /******* /**
/**////  /**    /**/**///**  /**      /**/**    *****/**////  /**
/**      /**    /**/**  //** //**     ** //**  ////**/**      /**
/********//******* /**   //** //*******   //******** /********/********
////////  ///////  //     //   ///////     ////////  //////// ////////
"""
# %%
''' de poner enero de 2004 he puesto 15 de enero de 2004, para todos los meses
'''
data0 = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/bq/set5/MEDIAS clorofila_histórico para Eugenio.xlsx',  dtype={  
"E2": "float64", "E6": "float64", "E13": "float64","E17": "float64","E23": "float64","E24": "float64",}, parse_dates= ['mes'], sheet_name='EUROGEL', skiprows=1, usecols= ['mes', 'E2', 'E6', 'E13', 'E17', 'E23','E24'], nrows=25) 
print(data0.head())
print(f'tiene una longitud de {len(data0)} filas')
# %%
data = data0.sort_values("mes").reset_index(drop=True)
# %%
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (data.mes - epoch) / pd.Timedelta(days=1)
# %%
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/set5/chl eugenio/'
ncfile = Dataset(f'{path}clorofila_eugenio_EUROGEL.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='clorofila eugenio pçagina EUROGEL'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'; ncfile.source = 'XXX'; ncfile.Conventions = 'CF-1.8'
# %%
max_length_param = len("E24")
estaciones = ['E2', 'E6', 'E13', 'E17','E23', 'E24']
estaciones_np = np.array(estaciones, dtype=f'S{max(len(s) for s in estaciones)}')

df_chl =  data.iloc[:, [1,2,3, 4,5,6]]
ncfile.createDimension('time', len(dias_desde_1970))
ncfile.createDimension('unit_char_len', max_length_param)
ncfile.createDimension('station_name', len(estaciones))
for dim in ncfile.dimensions.items():
    print(dim)

# %%
time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
time_var[:] = dias_desde_1970.values  # Se asigna directamente

parameter_var = ncfile.createVariable('station', 'S1', ('station_name', 'unit_char_len'))
parameter_var.long_name = 'station'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)
# %%
value_var = ncfile.createVariable('chlorophyll', np.float64, ('time','station_name'))
value_var.long_name = 'mass_concentration_of_chlorophyll_in_sea_water'
value_var.units= '¿mg L-1?'
value_var[:,:] = df_chl

lat_var = ncfile.createVariable('latitude', np.float64,('station_name') )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var[:] = [ 9999, 9999, 9999, 9999, 9999, 9999]

lon_var = ncfile.createVariable('longitude', np.float64,('station_name') )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var[:] = [ 9999, 9999, 9999 , 9999,  9999, 9999]

crs = ncfile.createVariable("crs", "i") 
crs.long_name=  "Reference coordinate system"
crs.grid_mapping_name = "latitude_longitude"
crs.epsg_code = "EPSG:4326"

# %%
ncfile.close()
# %% COMPROBACION
dataset = Dataset(f'{path}clorofila_eugenio_EUROGEL.nc', "r")
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
unit = dataset.variables["station"][:]    #  
value = dataset.variables["chlorophyll"][:] 
print(tiempo); print('-----------------')
print(unit); print('-----------------')
print(value); print('-----------------')
# %%
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()

# %% PLOT
fechas = np.array(fechas)  # convertir a numpy array para indexar
orden = np.argsort(fechas)  # obtener índices ordenados

# Reordenar fechas y nitrato
fechas_ordenadas = fechas[orden]
chl_ordenado = value[orden, :]

fig, axes = plt.subplots(nrows=len(estaciones), figsize=(10, 10), sharex=True)

for i in range(len(estaciones)):
    axes[i].plot(fechas_ordenadas, chl_ordenado[:, i], marker='o', label=estaciones[i])
    axes[i].set_ylabel('clorofila ')
    axes[i].set_title(f'clorofila_eugenio_EUROGEL- {estaciones[i]}')
    axes[i].grid(True)
    axes[i].legend()

axes[-1].set_xlabel('Fecha')
plt.tight_layout()
plt.show()

# %%
generar_txt(f'{path}clorofila_eugenio_EUROGEL.nc', f'{path}clorofila_eugenio_EUROGEL_display.txt')
# %%
"""
  ********     **     ********** ******** **       ** ********** ********
 **//////     ****   /////**/// /**///// /**      /**/////**/// /**/////
/**          **//**      /**    /**      /**      /**    /**    /**
/*********  **  //**     /**    /******* /**      /**    /**    /*******
////////** **********    /**    /**////  /**      /**    /**    /**////
       /**/**//////**    /**    /**      /**      /**    /**    /**
 ******** /**     /**    /**    /********/********/**    /**    /********
////////  //      //     //     //////// //////// //     //     ////////
"""
# %%
data0 = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/bq/set5/MEDIAS clorofila_histórico para Eugenio.xlsx',  dtype={  
"CHLA satelital corregida": "float64",}, parse_dates= ['FECHA'], sheet_name='SATELITE', skiprows=0, usecols= [ 'FECHA','CHLA satelital corregida'], ) 
print(data0.head())
print(f'tiene una longitud de {len(data0)} filas')
# %%
data = data0.sort_values("FECHA").reset_index(drop=True)
# %%
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (data.FECHA - epoch) / pd.Timedelta(days=1)
# %%
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/set5/chl eugenio/'
ncfile = Dataset(f'{path}clorofila_eugenio_SATELITE.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='clorofila eugenio pa  gina SATELITE'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'; ncfile.source = 'XXX'; ncfile.Conventions = 'CF-1.8'
# %%
ncfile.createDimension('time', len(dias_desde_1970))

for dim in ncfile.dimensions.items():
    print(dim)

# %%
time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
time_var[:] = dias_desde_1970.values  # Se asigna directamente

# %%
value_var = ncfile.createVariable('chlorophyll', np.float64, ('time',))
value_var.long_name = 'mass_concentration_of_chlorophyll_in_sea_water'
value_var.units= '¿mg L-1?'
value_var.comment = 'corrected satellited chla'
value_var[:] =data['CHLA satelital corregida']


lat_var = ncfile.createVariable('latitude', np.float64, )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var[:] = [ 9999]

lon_var = ncfile.createVariable('longitude', np.float64, )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var[:] = [ 9999]

crs = ncfile.createVariable("crs", "i") 
crs.long_name=  "Reference coordinate system"
crs.grid_mapping_name = "latitude_longitude"
crs.epsg_code = "EPSG:4326"
# %%
ncfile.close()
# %% COMPROBACION
dataset = Dataset(f'{path}clorofila_eugenio_SATELITE.nc', "r")
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
value = dataset.variables["chlorophyll"][:] 
print(tiempo); print('-----------------')
print(value); print('-----------------')
# %%
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()

# %% PLOT
fechas = np.array(fechas)  # convertir a numpy array para indexar
orden = np.argsort(fechas)  # obtener índices ordenados

# Reordenar fechas y nitrato
fechas_ordenadas = fechas[orden]
chl_ordenado = value[orden]

fig, axes = plt.subplots(figsize=(10, 10))

axes.plot(fechas_ordenadas, chl_ordenado, marker='o', label=estaciones[i])
axes.set_ylabel('clorofila ')
axes.set_title(f'clorofila_eugenio_SATELITE')
axes.grid(True)
axes.legend()

axes.set_xlabel('Fecha')
plt.tight_layout()
plt.show()

# %%
generar_txt(f'{path}clorofila_eugenio_SATELITE.nc', f'{path}clorofila_eugenio_SATELITE_display.txt')

# %%
