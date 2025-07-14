# %%
from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
from generar_txt import generar_txt
# %%
data  = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/bq/set3/Nutrient_Belich_2021.xlsx', parse_dates= ['Date'],sheet_name="data")
data = data.sort_values(by='Date').reset_index(drop=True)
print(data.head())

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
# path_viejo = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/set3-julia/nitrato/'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/nutrients/nitrate/BELICH_NUT'
ncfile = Dataset(f'{path}BELICH_NUT_NO3.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)
# C:\Users\lilia.flores\Nextcloud\Documents\Datos_MM_Art_2025\datasets_ncFormat\Biogeochemical\nutrients\nitrate\BELICH_NUT

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='BELICH_NUT_NO3'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'; ncfile.source = 'XXX'; ncfile.Conventions = 'CF-1.8'
ncfile.comment = 'A -1 value means LOWER THAN DETECTION LIMIT'
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
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
time_var[:] = dias_desde_1970.values  # Se asigna directamente

parameter_var = ncfile.createVariable('station_name', 'S1', ('station', 'unit_char_len'))
parameter_var.long_name = 'station'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)
# %%
value_var = ncfile.createVariable('nitrate', np.float64, ('time','station'))
value_var.long_name = 'mass concentration of nitrate in sea water'
value_var.units= '¿mg L-1?'
value_var[:,:] = df_nitratos

lat_var = ncfile.createVariable('latitude', np.float64,('station') )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var[:] = [ 37.792949, 37.698052, 37.667697]

lon_var = ncfile.createVariable('longitude', np.float64,('station') )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var[:] = [ -0.78331724, -0.78319145 ,-0.75235986]

crs = ncfile.createVariable("crs", "i")  # Dummy variable for coordinate reference system
crs.long_name=  "Reference coordinate system"
crs.grid_mapping_name = "latitude_longitude"
crs.epsg_code = "EPSG:4326"

# %%
ncfile.close()

# %% COMPROBACION
dataset = Dataset(f'{path}nutrientes_3_estaciones_NITRATO.nc', "r")
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
generar_txt(f'{path}BELICH_NUT_NO3.nc', f'{path}BELICH_NUT_NO3_display.txt')

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
# path_viejo= 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/set3-julia/nitrito/'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/nutrients/nitrite/BELICH_NUT'


ncfile = Dataset(f'{path}BELICH_NUT_NO2.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='BELICH_NUT_NO2'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'; ncfile.source = 'XXX'; ncfile.Conventions = 'CF-1.8'
ncfile.comments = 'A -1 value means LOWER THAN DETECTION LIMIT'

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
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
time_var[:] = dias_desde_1970.values  # Se asigna directamente

parameter_var = ncfile.createVariable('station_name', 'S1', ('station', 'unit_char_len'))
parameter_var.long_name = 'station'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)
# %%
value_var = ncfile.createVariable('nitrite', np.float64, ('time','station'))
value_var.long_name = 'mass concentration of nitrite in sea water'
value_var.units= '¿mg L-1?'
value_var[:,:] = df_nitritos

lat_var = ncfile.createVariable('latitude', np.float64,('station') )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var[:] = [ 37.792949, 37.698052, 37.667697]

lon_var = ncfile.createVariable('longitude', np.float64,('station') )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var[:] = [ -0.78331724, -0.78319145 ,-0.75235986]

crs = ncfile.createVariable("crs", "i")  # Dummy variable for coordinate reference system
crs.long_name=  "Reference coordinate system"
crs.grid_mapping_name = "latitude_longitude"
crs.epsg_code = "EPSG:4326"


# %%
ncfile.close()

# %% COMPROBACION
dataset = Dataset(f'{path}BELICH_NUT_NO2.nc', "r")
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
generar_txt(f'{path}BELICH_NUT_NO2.nc', f'{path}BELICH_NUT_NO2_display.txt')

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
# path_viejo= 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/set3-julia/fosfato/'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/nutrients/phosphate/BELICH_NUT'
ncfile = Dataset(f'{path}BELICH_NUT_PO4.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='BELICH_NUT_PO4'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'; ncfile.source = 'XXX'; ncfile.Conventions = 'CF-1.8'
ncfile.comment = 'A -1 value means LOWER THAN DETECTION LIMIT'

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
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
time_var[:] = dias_desde_1970.values  # Se asigna directamente

parameter_var = ncfile.createVariable('station_name', 'S1', ('station', 'unit_char_len'))
parameter_var.long_name = 'station'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)
# %%
value_var = ncfile.createVariable('phosphate', np.float64, ('time','station'))
value_var.long_name = 'mass concentration of phosphate in sea water'
value_var.units= '¿mg L-1?'
value_var[:,:] = df_fosfatos

lat_var = ncfile.createVariable('latitude', np.float64,('station') )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var[:] = [ 37.792949, 37.698052, 37.667697]

lon_var = ncfile.createVariable('longitude', np.float64,('station') )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var[:] = [ -0.78331724, -0.78319145 ,-0.75235986]

crs = ncfile.createVariable("crs", "i")
crs.long_name=  "Reference coordinate system"
crs.grid_mapping_name = "latitude_longitude"
crs.epsg_code = "EPSG:4326"
# %%
ncfile.close()

# %% COMPROBACION
dataset = Dataset(f'{path}BELICH_NUT_PO4.nc', "r")
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
generar_txt(f'{path}BELICH_NUT_PO4.nc', f'{path}BELICH_NUT_PO4_display.txt')
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
# path_viejo= 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/set3-julia/silicato/'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/nutrients/silicate/SiO4'
ncfile = Dataset(f'{path}BELICH_NUT_SiO4.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='BELICH_NUT_SiO4'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'; ncfile.source = 'XXX'; ncfile.Conventions = 'CF-1.8'
ncfile.comment = 'A -1 value means LOWER THAN DETECTION LIMIT'
# %%
df_silicatos =  data.iloc[:, [4,10,16]]
ncfile.createDimension('time', len(dias_desde_1970))
ncfile.createDimension('unit_char_len', max_length_param)
ncfile.createDimension('station', len(estaciones))
for dim in ncfile.dimensions.items():
    print(dim)

# %%
time_var = ncfile.createVariable('time', np.float64, ('time',))
time_var.units = "days since 1970-01-01 00:00:0"
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
time_var[:] = dias_desde_1970.values  # Se asigna directamente

parameter_var = ncfile.createVariable('station_name', 'S1', ('station', 'unit_char_len'))
parameter_var.long_name = 'station'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)
# %%
value_var = ncfile.createVariable('silicate', np.float64, ('time','station'))
value_var.long_name = 'mass concentration of silicate in sea water'
value_var.units= '¿mg L-1?'
value_var[:,:] = df_silicatos

lat_var = ncfile.createVariable('latitude', np.float64,('station') )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var[:] = [ 37.792949, 37.698052, 37.667697]

lon_var = ncfile.createVariable('longitude', np.float64,('station') )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var[:] = [ -0.78331724, -0.78319145 ,-0.75235986]

crs = ncfile.createVariable("crs", "i")
crs.long_name = "Reference coordinate system"
crs.grid_mapping_name = "latitude_longitude"
crs.epsg_code = "EPSG:4326"
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
generar_txt(f'{path}BELICH_NUT_SiO4.nc', f'{path}BELICH_NUT_SiO4_display.txt')
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
# path_viejo = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/set3-julia/amonio/'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/nutrients/ammonium/NH4'
ncfile = Dataset(f'{path}BELICH_NUT_NH4.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='BELICH_NUT_NH4'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'; ncfile.source = 'XXX'; ncfile.Conventions = 'CF-1.8'
ncfile.comment = 'A -1 value means LOWER THAN DETECTION LIMIT'

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
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
time_var[:] = dias_desde_1970.values  # Se asigna directamente

parameter_var = ncfile.createVariable('station_name', 'S1', ('station', 'unit_char_len'))
parameter_var.long_name = 'station'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)
# %%
value_var = ncfile.createVariable('ammonium', np.float64, ('time','station'))
value_var.long_name = 'mass_concentration_of_ammonium_in_sea_water'
value_var.units= '¿mg L-1?'
value_var[:,:] = df_amonio

lat_var = ncfile.createVariable('latitude', np.float64,('station') )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var[:] = [ 37.792949, 37.698052, 37.667697]

lon_var = ncfile.createVariable('longitude', np.float64,('station') )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var[:] = [ -0.78331724, -0.78319145 ,-0.75235986]

crs = ncfile.createVariable("crs", "i")
crs.long_name=  "Reference coordinate system"
crs.grid_mapping_name = "latitude_longitude"
crs.epsg_code = "EPSG:4326"
# %%
ncfile.close()

# %% COMPROBACION
dataset = Dataset(f'{path}BELICH_NUT_NH4.nc', "r")
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
generar_txt(f'{path}BELICH_NUT_NH4.nc', f'{path}BELICH_NUT_NH4_display.txt')

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
# path_viejo = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/set3-julia/din/'
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/nutrients/din/BELICH_NUT_DIN'
ncfile = Dataset(f'{path}BELICH_NUT_DIN.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title='BELICH_NUT_DIN'
ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'XXXX'; ncfile.source = 'XXX'; ncfile.Conventions = 'CF-1.8'
ncfile.comment = 'A -1 value means LOWER THAN DETECTION LIMIT'

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
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
time_var[:] = dias_desde_1970.values  # Se asigna directamente

parameter_var = ncfile.createVariable('station_name', 'S1', ('station', 'unit_char_len'))
parameter_var.long_name = 'station'
parameter_var._Encoding = 'ascii'
parameter_var[:,:] = stringtochar(estaciones_np)

# %%
value_var = ncfile.createVariable('din', np.float64, ('time','station'))
value_var.long_name = 'mass_concentration_of_dissolved_inorganic_nitrogen_in_sea_water'
value_var.units= '¿mg -L?'
value_var[:,:] = df_din
value_var.comment = 'DIN is the sum of nitrate, nitrite and ammonium concentrations in seawater'

lat_var = ncfile.createVariable('latitude', np.float64,('station') )
lat_var.units = 'degrees_north'
lat_var.standard_name = "latitude"
lat_var[:] = [ 37.792949, 37.698052, 37.667697]

lon_var = ncfile.createVariable('longitude', np.float64,('station') )
lon_var.units = 'degrees_east'
lon_var.standard_name = "longitude"
lon_var[:] = [ -0.78331724, -0.78319145 ,-0.75235986]

crs = ncfile.createVariable("crs", "i")
crs.long_name=  "Reference coordinate system"
crs.grid_mapping_name = "latitude_longitude"
crs.epsg_code = "EPSG:4326"
# %%
ncfile.close()

# %% COMPROBACION
dataset = Dataset(f'{path}BELICH_NUT_DIN.nc', "r")
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
generar_txt(f'{path}BELICH_NUT_DIN.nc', f'{path}BELICH_NUT_DIN.txt')
# %%
