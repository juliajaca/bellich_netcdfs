# %%
from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
from generar_txt import generar_txt
# %%
data  = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/bq/set3/Nutrient_Belich_2021.xlsx', parse_dates= ['Date'],sheet_name="data")
data = data.sort_values(by='Date').reset_index(drop=True)
print(data.head())
data = data.replace(-1, np.nan)
# %%
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (data.Date - epoch) / pd.Timedelta(days=1)
# %%
"""
          **   **                      **
         //   /**                     /**
 *******  ** ****** ******  ******   ******  ******
//**///**/**///**/ //**//* //////** ///**/  **////**
 /**  /**/**  /**   /** /   *******   /**  /**   /**
 /**  /**/**  /**   /**    **////**   /**  /**   /**
 ***  /**/**  //** /***   //********  //** //******
///   // //    //  ///     ////////    //   //////
"""

nombre_fichero = 'BELICH_NUT_NO3'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/nutrients/nitrate/BELICH_NUT/'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)


# %% CREAR ATRIBUTOS GLOBALES
ncfile.title= nombre_fichero
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'BELICH_NUT'
ncfile.project = 'BELICH'; ncfile.source = 'In situ data collection'; ncfile.Conventions = 'CF-1.8'
ncfile.comment = 'Values of -1, indicating concentrations below the detection limit, have been replaced with NaN'

# %%
# crear dimensiones
max_length_param = len("station X")
estaciones = ['Station A', 'Station B', 'Station C']
estaciones_np = np.array(estaciones, dtype=f'S{max(len(s) for s in estaciones)}')

df_nitratos =  data.iloc[:, [1,7,13]]
ncfile.createDimension('time', len(dias_desde_1970))
ncfile.createDimension('unit_char_len', max_length_param)
ncfile.createDimension('station', len(estaciones))
for dim in ncfile.dimensions.items():
    print(dim)

# %%
time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values  # Se asigna directamente

lat_var = ncfile.createVariable('latitude', np.float64,('station') )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var.grid_mapping = "crs"
lat_var[:] = [ 37.792949, 37.698052, 37.667697]

lon_var = ncfile.createVariable('longitude', np.float64,('station') )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] = [ -0.78331724, -0.78319145 ,-0.75235986]

# %%
print(df_nitratos)
valores_con_nan = ( df_nitratos *1000) / 62.0049 #paso de mg L-1 a umol L-1 
valores_con_nan[np.isnan(valores_con_nan)] = -9999
print(valores_con_nan)
value_var = ncfile.createVariable('nitrate', np.float64, ('time','station'))
value_var.units= 'umol L-1'
value_var.standard_name = "mole_concentration_of_nitrate_in_sea_water"
value_var.long_name = 'Nitrate concentration in sea water'
value_var.cell_methods = "time: mean"
value_var.missing_value = -9999
value_var.grid_mapping = "crs"
value_var[:,:] = valores_con_nan
value_var.comment  = "converted from original data assumed to be in mg/L to µmol/L "

parameter_var = ncfile.createVariable('station_name', 'S1', ('station', 'unit_char_len'))
parameter_var.long_name = 'station'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)

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
unit = dataset.variables["station_name"][:]    #  
value = dataset.variables["nitrate"][:] 
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
nitrato_ordenado = value[orden, :]

fig, axes = plt.subplots(nrows=3, figsize=(10, 10), sharex=True)

for i in range(3):
    axes[i].plot(fechas_ordenadas, nitrato_ordenado[:, i], marker='o', label=estaciones[i])
    axes[i].set_ylabel('Nitrato ')
    axes[i].set_title(f'Serie temporal NITRATO- {estaciones[i]}')
    axes[i].grid(True)
    axes[i].legend()

axes[-1].set_xlabel('Fecha')
plt.tight_layout()
plt.show()

# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}.txt')
# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Biogeochemical/nutrients/nitrate/BELICH_NUT/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')
"""
.##....##.####.########.########..####.########..#######.
.###...##..##.....##....##.....##..##.....##....##.....##
.####..##..##.....##....##.....##..##.....##....##.....##
.##.##.##..##.....##....########...##.....##....##.....##
.##..####..##.....##....##...##....##.....##....##.....##
.##...###..##.....##....##....##...##.....##....##.....##
.##....##.####....##....##.....##.####....##.....#######.
"""
# %%

path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/nutrients/nitrite/BELICH_NUT/'
nombre_fichero = 'BELICH_NUT_NO2'

ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title= nombre_fichero
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'BELICH_NUT'
ncfile.project = 'BELICH'; ncfile.source = 'In situ data collection'; ncfile.Conventions = 'CF-1.8'
ncfile.comment = 'Values of -1, indicating concentrations below the detection limit, have been replaced with NaN'

# %%
df_nitritos =  data.iloc[:, [2,8,14]]
ncfile.createDimension('time', len(dias_desde_1970))
ncfile.createDimension('unit_char_len', max_length_param)
ncfile.createDimension('station', len(estaciones))
for dim in ncfile.dimensions.items():
    print(dim)

# %%
time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values  # Se asigna directamente

lat_var = ncfile.createVariable('latitude', np.float64,('station') )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var.grid_mapping = "crs"
lat_var[:] = [ 37.792949, 37.698052, 37.667697]

lon_var = ncfile.createVariable('longitude', np.float64,('station') )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] = [ -0.78331724, -0.78319145 ,-0.75235986]

# %%
print(df_nitritos)
valores_con_nan = ( df_nitritos *1000) / 46.01 #paso de mg L-1 a umol L-1 
valores_con_nan[np.isnan(valores_con_nan)] = -9999
print(valores_con_nan)
value_var = ncfile.createVariable('nitrite', np.float64, ('time','station'))
value_var.units= 'umol L-1'
value_var.standard_name = "mole_concentration_of_nitrite_in_sea_water"
value_var.long_name = 'Nitrite concentration in sea water'
value_var.cell_methods = "time: mean"
value_var.missing_value = -9999
value_var.grid_mapping = "crs"
value_var.comment  = "converted from original data assumed to be in mg/L to µmol/L "
value_var[:,:] = valores_con_nan

parameter_var = ncfile.createVariable('station_name', 'S1', ('station', 'unit_char_len'))
parameter_var.long_name = 'station'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)

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
unit = dataset.variables["station_name"][:]    #  
value = dataset.variables["nitrite"][:] 
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
nitrito_ordenado = value[orden, :]

fig, axes = plt.subplots(nrows=3, figsize=(10, 10), sharex=True)

for i in range(3):
    axes[i].plot(fechas_ordenadas, nitrito_ordenado[:, i], marker='o', label=estaciones[i])
    axes[i].set_ylabel('Nitrito ')
    axes[i].set_title(f'Serie temporal NITRITO- {estaciones[i]}')
    axes[i].grid(True)
    axes[i].legend()

axes[-1].set_xlabel('Fecha')
plt.tight_layout()
plt.show()

# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}.txt')

# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Biogeochemical/nutrients/nitrite/BELICH_NUT/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')
# %%
"""
 ********   *******    ******** ********     **     **********   *******
/**/////   **/////**  **////// /**/////     ****   /////**///   **/////**
/**       **     //**/**       /**         **//**      /**     **     //**
/******* /**      /**/*********/*******   **  //**     /**    /**      /**
/**////  /**      /**////////**/**////   **********    /**    /**      /**
/**      //**     **        /**/**      /**//////**    /**    //**     **
/**       //*******   ******** /**      /**     /**    /**     //*******
//         ///////   ////////  //       //      //     //       ///////
"""
# %%
nombre_fichero ='BELICH_NUT_PO4'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/nutrients/phosphate/BELICH_NUT'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title= nombre_fichero
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'BELICH_NUT'
ncfile.project = 'BELICH'; ncfile.source = 'In situ data collection'; ncfile.Conventions = 'CF-1.8'
ncfile.comment = 'Values of -1, indicating concentrations below the detection limit, have been replaced with NaN'

# %%
df_fosfatos =  data.iloc[:, [3,9,15]]
ncfile.createDimension('time', len(dias_desde_1970))
ncfile.createDimension('unit_char_len', max_length_param)
ncfile.createDimension('station', len(estaciones))
for dim in ncfile.dimensions.items():
    print(dim)

# %%
time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values  # Se asigna directamente

# %%
lat_var = ncfile.createVariable('latitude', np.float64,('station') )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var.grid_mapping = "crs"
lat_var[:] = [ 37.792949, 37.698052, 37.667697]

lon_var = ncfile.createVariable('longitude', np.float64,('station') )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] = [ -0.78331724, -0.78319145 ,-0.75235986]

print(df_fosfatos)
valores_con_nan = ( df_fosfatos *1000) / 94.97 #paso de mg L-1 a umol L-1 
valores_con_nan[np.isnan(valores_con_nan)] = -9999
print(valores_con_nan)

value_var = ncfile.createVariable('phosphate', np.float64, ('time','station'))
value_var.units= 'umol L-1'
value_var.standard_name= 'mole_concentration_of_phosphate_in_sea_water'
value_var.long_name = 'Phosphate concentration in sea water'
value_var.cell_methods = "time: mean"
value_var.missing_value = -9999
value_var.grid_mapping = "crs"
value_var.comment  = "Converted from original data assumed to be in mg/L to µmol/L "
value_var[:,:] = valores_con_nan

parameter_var = ncfile.createVariable('station_name', 'S1', ('station', 'unit_char_len'))
parameter_var.long_name = 'station'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)

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
unit = dataset.variables["station_name"][:]    #  
value = dataset.variables["phosphate"][:] 
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
nitrito_ordenado = value[orden, :]

fig, axes = plt.subplots(nrows=3, figsize=(10, 10), sharex=True)

for i in range(3):
    axes[i].plot(fechas_ordenadas, nitrito_ordenado[:, i], marker='o', label=estaciones[i])
    axes[i].set_ylabel('Fosfato ')
    axes[i].set_title(f'Serie temporal FOSFATO- {estaciones[i]}')
    axes[i].grid(True)
    axes[i].legend()

axes[-1].set_xlabel('Fecha')
plt.tight_layout()
plt.show()

# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}.txt')
# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Biogeochemical/nutrients/phosphate/BELICH_NUT/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')
# %%
"""
  ******** ** **       **   ******      **     **********   *******
 **////// /**/**      /**  **////**    ****   /////**///   **/////**
/**       /**/**      /** **    //    **//**      /**     **     //**
/*********/**/**      /**/**         **  //**     /**    /**      /**
////////**/**/**      /**/**        **********    /**    /**      /**
       /**/**/**      /**//**    **/**//////**    /**    //**     **
 ******** /**/********/** //****** /**     /**    /**     //*******
////////  // //////// //   //////  //      //     //       ///////
"""
# %%
nombre_fichero ='BELICH_NUT_SiO4'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/nutrients/silicate/'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title= nombre_fichero
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'BELICH_NUT'
ncfile.project = 'BELICH'; ncfile.source = 'In situ data collection'; ncfile.Conventions = 'CF-1.8'
ncfile.comment = 'Values of -1, indicating concentrations below the detection limit, have been replaced with NaN'
# %%
df_silicatos =  data.iloc[:, [4,10,16]]
print(df_silicatos)
valores_con_nan = ( df_silicatos *1000) / 96.06 #paso de mg L-1 a umol L-1 
valores_con_nan[np.isnan(valores_con_nan)] = -9999
print(valores_con_nan)

ncfile.createDimension('time', len(dias_desde_1970))
ncfile.createDimension('unit_char_len', max_length_param)
ncfile.createDimension('station', len(estaciones))
for dim in ncfile.dimensions.items():
    print(dim)

time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values  # Se asigna directamente


# %%
lat_var = ncfile.createVariable('latitude', np.float64,('station') )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var.grid_mapping = "crs"
lat_var[:] = [ 37.792949, 37.698052, 37.667697]

lon_var = ncfile.createVariable('longitude', np.float64,('station') )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] = [ -0.78331724, -0.78319145 ,-0.75235986]

value_var = ncfile.createVariable('silicate', np.float64, ('time','station'))
value_var.units= 'umol L-1'
value_var.standard_name = "mole_concentration_of_silicate_in_sea_water"
value_var.long_name = 'Silicate concentration in sea water'
value_var.cell_methods = "time: mean"
value_var.missing_value = -9999
value_var.grid_mapping = "crs"
value_var.comment  = "Converted from original data assumed to be in mg/L to µmol/L "
value_var[:,:] = valores_con_nan

parameter_var = ncfile.createVariable('station_name', 'S1', ('station', 'unit_char_len'))
parameter_var.long_name = 'station'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)
# %%
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

# %% COMPROBACION
dataset = Dataset(f'{path}BELICH_NUT_SiO4.nc', "r")
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
unit = dataset.variables["station_name"][:]    #  
value = dataset.variables["silicate"][:] 
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
silicato_ordenado = value[orden, :]

fig, axes = plt.subplots(nrows=3, figsize=(10, 10), sharex=True)

for i in range(3):
    axes[i].plot(fechas_ordenadas, silicato_ordenado[:, i], marker='o', label=estaciones[i])
    axes[i].set_ylabel('Silicato')
    axes[i].set_title(f'Serie temporal SILICATO- {estaciones[i]}')
    axes[i].grid(True)
    axes[i].legend()

axes[-1].set_xlabel('Fecha')
plt.tight_layout()
plt.show()

# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}.txt')

# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Biogeochemical/nutrients/silicate/BELICH_NUT/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')
# %%
"""
     **     ****     ****   *******   ****     ** **   *******
    ****   /**/**   **/**  **/////** /**/**   /**/**  **/////**
   **//**  /**//** ** /** **     //**/**//**  /**/** **     //**
  **  //** /** //***  /**/**      /**/** //** /**/**/**      /**
 **********/**  //*   /**/**      /**/**  //**/**/**/**      /**
/**//////**/**   /    /**//**     ** /**   //****/**//**     **
/**     /**/**        /** //*******  /**    //***/** //*******
//      // //         //   ///////   //      /// //   ///////
"""
nombre_fichero='BELICH_NUT_NH4'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/nutrients/ammonium/NH4'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title= nombre_fichero
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'BELICH_NUT'
ncfile.project = 'BELICH'; ncfile.source = 'In situ data collection'; ncfile.Conventions = 'CF-1.8'
ncfile.comment = 'Values of -1, indicating concentrations below the detection limit, have been replaced with NaN'

# %%
df_amonio =  data.iloc[:, [5,11,17]]
ncfile.createDimension('time', len(dias_desde_1970))
ncfile.createDimension('unit_char_len', max_length_param)
ncfile.createDimension('station', len(estaciones))
for dim in ncfile.dimensions.items():
    print(dim)

# %%
time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values  # Se asigna directamente

lat_var = ncfile.createVariable('latitude', np.float64,('station') )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var.grid_mapping = "crs"
lat_var[:] = [ 37.792949, 37.698052, 37.667697]

lon_var = ncfile.createVariable('longitude', np.float64,('station') )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] = [ -0.78331724, -0.78319145 ,-0.75235986]
# %%
print(df_amonio)
valores_con_nan = ( df_amonio *1000) / 18.04 #paso de mg L-1 a umol L-1 
valores_con_nan[np.isnan(valores_con_nan)] = -9999
print(valores_con_nan)

value_var = ncfile.createVariable('ammonium', np.float64, ('time','station'))
value_var.units= 'umol L-1'
value_var.standard_name = "mole_concentration_of_ammonium_in_sea_water"
value_var.long_name = 'Ammonium concentration in sea water'
value_var.cell_methods = "time: mean"
value_var.missing_value = -9999
value_var.grid_mapping = "crs"
value_var.comment  = "Converted from original data assumed to be in mg/L to µmol/L "
value_var[:,:] = valores_con_nan

parameter_var = ncfile.createVariable('station_name', 'S1', ('station', 'unit_char_len'))
parameter_var.long_name = 'station'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)

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
unit = dataset.variables["station_name"][:]    #  
value = dataset.variables["ammonium"][:] 
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
amonio_ordenado = value[orden, :]

fig, axes = plt.subplots(nrows=3, figsize=(10, 10), sharex=True)

for i in range(3):
    axes[i].plot(fechas_ordenadas, amonio_ordenado[:, i], marker='o', label=estaciones[i])
    axes[i].set_ylabel('Amonio ')
    axes[i].set_title(f'Serie temporal AMONIO- {estaciones[i]}')
    axes[i].grid(True)
    axes[i].legend()

axes[-1].set_xlabel('Fecha')
plt.tight_layout()
plt.show()

# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}.txt')
# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Biogeochemical/nutrients/ammonium/BELICH_NUT/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')

# %% 
"""
 *******   ** ****     **
/**////** /**/**/**   /**
/**    /**/**/**//**  /**
/**    /**/**/** //** /**
/**    /**/**/**  //**/**
/**    ** /**/**   //****
/*******  /**/**    //***
///////   // //      ///
"""
# %%
nombre_fichero ='BELICH_NUT_DIN'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/nutrients/din/N'
ncfile = Dataset(f'{path}{nombre_fichero}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title= nombre_fichero
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'BELICH_NUT'
ncfile.project = 'BELICH'; ncfile.source = 'In situ data collection'; ncfile.Conventions = 'CF-1.8'
ncfile.comment = 'Values of -1, indicating concentrations below the detection limit, have been replaced with NaN'

# %%
df_din = data.iloc[:, [6,12,18]]
ncfile.createDimension('time', len(dias_desde_1970))
ncfile.createDimension('unit_char_len', max_length_param)
ncfile.createDimension('station', len(estaciones))
for dim in ncfile.dimensions.items():
    print(dim)

# %%
time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.calendar = 'gregorian'
time_var.standard_name = "time"
time_var[:] = dias_desde_1970.values  # Se asigna directamente

lat_var = ncfile.createVariable('latitude', np.float64,('station') )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var.grid_mapping = "crs"
lat_var[:] = [ 37.792949, 37.698052, 37.667697]

lon_var = ncfile.createVariable('longitude', np.float64,('station') )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var.grid_mapping = "crs"
lon_var[:] = [ -0.78331724, -0.78319145 ,-0.75235986]

print(df_amonio)
valores_con_nan = ( df_amonio *1000) / 14.01 #paso de mg L-1 a umol L-1 
valores_con_nan[np.isnan(valores_con_nan)] = -9999
print(valores_con_nan)

value_var = ncfile.createVariable('din', np.float64, ('time','station'))
value_var.units= 'umol L-1'
value_var.standard_name = "mole_concentration_of_dissolved_inorganic_nitrogen_in_sea_water"
value_var.long_name = 'Dissolved inorganic nitrogen concentration in sea water'
value_var.cell_methods = "time: mean"
value_var.missing_value = -9999
value_var.grid_mapping = "crs"
value_var.comment  = "Converted from original data assumed to be in mg/L to µmol/L "
value_var[:,:] = valores_con_nan

parameter_var = ncfile.createVariable('station_name', 'S1', ('station', 'unit_char_len'))
parameter_var.long_name = 'station'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)

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
unit = dataset.variables["station_name"][:]    #  
value = dataset.variables["din"][:] 
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
din_ordenado = value[orden, :]

fig, axes = plt.subplots(nrows=3, figsize=(10, 10), sharex=True)

for i in range(3):
    axes[i].plot(fechas_ordenadas,din_ordenado[:, i], marker='o', label=estaciones[i])
    axes[i].set_ylabel('din')
    axes[i].set_title(f'Serie temporal DIN- {estaciones[i]}')
    axes[i].grid(True)
    axes[i].legend()

axes[-1].set_xlabel('Fecha')
plt.tight_layout()
plt.show()

# %%
generar_txt(f'{path}{nombre_fichero}.nc', f'{path}{nombre_fichero}.txt')
# %%
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Biogeochemical/nutrients/din/BELICH_NUT/'
shutil.copy(f'{path}{nombre_fichero}.nc',f'{ruta_destino}{nombre_fichero}.nc')
# %%
