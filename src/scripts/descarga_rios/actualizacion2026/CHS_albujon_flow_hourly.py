# %%
from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import shutil
import sys
import pandas as pd
import matplotlib.pyplot as plt
sys.path.append("../../")
from generar_txt import generar_txt

# %%
data = pd.read_csv('C:/Users/Julia/Documents/VSCODE_BELLICH/src/scripts/descarga_rios/actualizacion2026/datos/CHS/CHS_Albujon_flow_hourly.csv', parse_dates= ['Date_time'], )

data["Discharge"] = pd.to_numeric(data["Discharge"].replace('#¡VALOR!', pd.NA), errors="coerce")
data["Date_time"] = pd.to_datetime(data["Date_time"], dayfirst=True)

nombre_fichero= 'CHS_Albujon_FLOW_HOURLY'
# %%
#  conversión de fechas : Fecha de referencia (1 de enero de 1970)
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (data.Date_time - epoch) / pd.Timedelta(days=1)

# %%
path = f'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Runoff/CHS_Albujon/{nombre_fichero}.nc'

ncfile = Dataset(path, mode='w', format='NETCDF3_CLASSIC')
print(ncfile)

# %% CREAR ATRIBUTOS GLOBALES
ncfile.title=f'{nombre_fichero}'
ncfile.institution="Confederación Hidrográfica del Segura (CHS)"
ncfile.domain= 'Mar menor coastal lagoon, Spain'
ncfile.dataset_id = 'CHS'
ncfile.project = 'Not associated with a specific project'
ncfile.source = 'In situ data collection'
ncfile.Conventions = "CF-1.8"
ncfile.comment = ('''
    Continuous sensor recording. Water level measured with radar time-of-flight sensor Endress+Hauser Micropilot FMR20. Flow is the average of one measurement per second, and highest resolution data are offered as 5-minute means. Note that the value is not necessarily the flow discharging in the Mar Menor, as the gauge is located 200 m from the actual outlet into the lagoon and is not accounting the water diverted to El Albujón pump station to conduct the water to the desalination plant of El Mojón (16 km NE). In general, most of the days referred to in these data the pumping was not operating. This effect is negligible in high flows. Time stamps correspond to civil time, the way of adjustment during the times civil time changes from summer to winter standard and vice versa is unknown.
'''
)

# Crear dimensiones
ncfile.createDimension('time', len(data))

for dim in ncfile.dimensions.items():
    print(dim)
# %%
date_var = ncfile.createVariable('time', np.float64, ('time'))
date_var.units= "days since 1970-01-01 00:00:00"
date_var.calendar = 'gregorian'
date_var.standard_name = "time"
date_var[:] = dias_desde_1970.values

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

flow_var = ncfile.createVariable('water_volume_transport', np.float64, ('time',))
flow_var.units = 'm3 s-1'
flow_var.long_name = 'Hourly Mean River Discharge'
valores_con_nan = data["Flow_mean "].values
valores_con_nan[np.isnan(valores_con_nan)] = -9999 
flow_var.cell_methods= "time: mean"
flow_var[:] =  valores_con_nan
flow_var.missing_value = -9999 
flow_var.grid_mapping = "crs"

valores_con_nan= data["Discharge"].values
valores_con_nan[np.isnan(valores_con_nan)] = -9999 
disch_var = ncfile.createVariable('accumulated_water_volume_transport', np.float64, ('time',))
disch_var.units = 'm3'
disch_var.long_name= 'Hourly Accumulated River Discharge'
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
generar_txt(f'{path}', f'{path}_display.txt')

# %%
#  COPIAR
ruta_destino = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Runoff/'
shutil.copy(f'{path}',f'{ruta_destino}{nombre_fichero}.nc')


# %% COMPROBACION
dataset = Dataset(f'{path}', "r")
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
disch = dataset.variables['accumulated_water_volume_transport'][:] 
flow = dataset.variables['water_volume_transport'][:]


print(tiempo); print('-----------------')
print(disch); print('-----------------')
print(flow[0:10]); print('-----------------')

fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
dataset.close()

# %% PLOT
fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")

orden = np.argsort(fechas)
fechas = fechas.values[orden]
flow = flow[orden]
disch = disch[orden]

fig, axs = plt.subplots(nrows=2, ncols=1,figsize=(18, 10),sharex=True)
axs[0].plot(fechas, flow,color="tab:blue",linewidth=1)
axs[0].set_ylabel("Flow (m³ s⁻¹)")
axs[0].set_title("Hourly Mean River Flow")
axs[0].grid(True)

axs[1].plot( fechas,disch,color="tab:orange",linewidth=1)
axs[1].set_ylabel("Hourly Discharge (m³)")
axs[1].set_title("Accumulated Daily River Discharge")
axs[1].grid(True)

axs[1].set_xlabel("Fecha")

plt.tight_layout()
plt.xticks(rotation=45)

plt.savefig(
    f"C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Runoff/CHS_Albujon/hourly_flow_discharge_timeseries.png",
    dpi=300, bbox_inches="tight")

plt.show()

# %%