# %%
from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
from generar_txt import generar_txt
# %%
data0 = pd.read_csv('datos_satelite_limpios_2.csv',  dtype={ "clorofila": "float64", "latitud": "float64","longitud": "float64",}, parse_dates= ['fecha'])
data = data0.sort_values(by=["fecha"]).reset_index(drop=True)
# %%
fechas_unicas = np.sort(data["fecha"].unique())
fechas_unicas_ts = pd.to_datetime(fechas_unicas)
latitudes = np.sort(data["latitud"].unique())
longitudes = np.sort(data["longitud"].unique())

epoch = pd.Timestamp('1970-01-01')
# %%
# Calcular la diferencia en días
dias_desde_1970 = (fechas_unicas_ts - epoch) / pd.Timedelta(days=1)

# %%
nombre_archivo = 'IEO_SAT_CHL'
# %%
ncfile = Dataset(f'{nombre_archivo}.nc', mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
# ncfile.title=f'{nombre_archivo}'
ncfile.title=f'sentinel_3_clorofila_bela'

ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
ncfile.domain= 'Mar menor coastal lagoon'
ncfile.project = 'IEO_SAT'
ncfile.source = 'Copernicus'
ncfile.conventions = 'CF-1.8'


# Crear dimensiones
ncfile.createDimension("time", len(fechas_unicas)) 
ncfile.createDimension("latitud", len(latitudes)) 
ncfile.createDimension("longitud", len(longitudes))

for dim in ncfile.dimensions.items():
    print(dim)
# %%
time_var = ncfile.createVariable('time', np.float64, ('time'))
time_var.units= "days since 1970-01-01 00:00:0"
time_var.standard_name = "time"
time_var.calendar = 'gregorian'
time_var[:] = dias_desde_1970.values

lat_var = ncfile.createVariable('latitud', np.float64, ('latitud',))
lat_var.units = 'degrees north'
lat_var.standard_name = 'latitud'
lat_var[:] = latitudes

lon_var = ncfile.createVariable('longitud', np.float64, ('longitud',))
lon_var.units = 'degrees east'
lon_var.standard_name = 'longitud'
lon_var[:] = longitudes

chl_var = ncfile.createVariable('chlorophyll', np.float64, ("time", "latitud", "longitud"), fill_value=None)
chl_var.units = 'mg m-3'
chl_var.standard_name = 'chlorophyll_concentration_in_sea_water'
chl_var.comment = 'BELA chlorophyll concentration algorithm used to estimate CHL.'
chl_var.cell_methods = (
    '1. Spatial mesh definition:\n'
    'Calculates the average spatial resolution (distance between consecutive points).\n'
    'Generate a regular grid (lat/lon).\n'
    'Applies a geographic mask (based on bathymetry and coastline) to delimit the area of interest.\n'
    '2. Fill in nans for each date:\n'
    'Clears duplicates by coordinates.\n'
    'Median distances in latitude and longitude have been calculated and if there is a jump between values of more than 1.5 * median, a lat lon point with nan value was created at a median distance.\n'
    '3. Mesh creation to interpolate:\n'    
    'Interpolates the chlorophyll values on the generated mesh using the nearest method.\n'
    'Uses a cKDTree to calculate the minimum distance from each mesh point to the real data.\n'
    'Applies a distance mask to eliminate unreliable interpolations (based on the median of real distances * 0.75).'
)

# %%
# Escribir datos en el NetCDF (fila por fila)
for i, fecha in enumerate(fechas_unicas_ts):
    sub_df = data[data["fecha"] == fecha] # Filtrar datos de la fecha actual

    # Crear una matriz vacía para la clorofila
    chl_matrix = np.full((len(latitudes), len(longitudes)), np.nan, dtype=np.float64)
    
    # Llenar la matriz con los valores de clorofila
    for j, (lat, lon, chl) in enumerate(zip(sub_df["latitud"], sub_df["longitud"], sub_df["clorofila"])):
        lat_idx = np.where(latitudes == lat)[0][0]
        lon_idx = np.where(longitudes == lon)[0][0]
        chl_matrix[lat_idx, lon_idx] = chl  # Asignar valor de clorofila
    
    chl_var[i, :, :] = chl_matrix # Guardar la matriz en el NetCDF sin cargar todo en RAM

ncfile.close()

# %%
dataset = Dataset(f'{nombre_archivo}.nc', "r")
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
lat = dataset.variables["latitud"][:] 
lon = dataset.variables["longitud"][:] 
value = dataset.variables["chlorophyll"][:,:,:] 
print(tiempo); print('-----------------')
print(lat); print('-----------------')
print(lon); print('-----------------')
print(value); print('-----------------')
# %%
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()
# %% PLOT
lon2d, lat2d = np.meshgrid(lon, lat)
# Aplanar para scatter
lon_flat = lon2d.flatten()
lat_flat = lat2d.flatten()
cloro_flat = value[50,:,:].flatten()

fig, axes = plt.subplots(1, 1, figsize=(14, 6), sharex=True, sharey=True)
sc1=axes.scatter(lon_flat, lat_flat, c=cloro_flat,)
axes.set_xlabel("Longitud")
axes.set_ylabel("Latitud")
axes.grid(True)
axes.legend()
fig.colorbar(sc1, ax=axes, label='Clorofila')

plt.tight_layout()
plt.show()

# %%
generar_txt(f'{nombre_archivo}.nc', f'{nombre_archivo}_display.txt')
# %%
