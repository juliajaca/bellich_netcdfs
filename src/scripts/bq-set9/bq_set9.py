# %%
from netCDF4 import Dataset, stringtochar
import numpy as np
import pandas as pd
import shutil
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
from generar_txt import generar_txt
# %%
n_dataset ='BELICH_BIOGQ_V2'
dict = {
    'clorofila (espectofotometria)':{'nombre_fichero':f'{n_dataset}_CHL',
                                    'v_name': 'chlorophyll',
                                     'standard_name': 'mass_concentration_of_chlorophyll_in_sea_water',
                                     'long_name': 'Chlorophyll-a Concentration in Sea Water',
                                     'comment':'Chlorophyll-a concentration determined via spectrophotometry following -standard/method name ???-'
                                     'Original data assumed to be in mg m-3 ',
                                     'units': 'ug L-1',
                                'path':'chlorophyll'},

    'nitrato':{'v_name': 'nitrate',
               'nombre_fichero':f'{n_dataset}_NO3',
               'standard_name':'mole_concentration_of_nitrate_in_sea_water',
                'long_name':'Nitrate concentration in sea water',
                  'units': 'umol L-1',
                'comment': 'Converted from original data assumed to be in mg/L to µmol/L',
                'path': 'nutrients/nitrate'},

    'nitrito':{'v_name': 'nitrite',
               'nombre_fichero':f'{n_dataset}_NO2',
               'standard_name':'mole_concentration_of_nitrite_in_sea_water',
            'long_name': 'Nitrite concentration in sea water',
              'units': 'umol L-1',
            'comment': 'Converted from original data assumed to be in mg/L to µmol/L',
             'path': 'nutrients/nitrite'},

    'fosfato':{'v_name': 'phosphate',
               'nombre_fichero':f'{n_dataset}_PO4',
               'standard_name': 'mole_concentration_of_phosphate_in_sea_water',
            'long_name': 'Phosphate concentration in sea water',
            'units': 'umol L-1',
            'comment': 'Converted from original data assumed to be in mg/L to µmol/L',
             'path': 'nutrients/phosphate'},

    'silicato':{'v_name': 'silicate',
               'nombre_fichero':f'{n_dataset}_SiO4',
                'standard_name':'mole_concentration_of_silicate_in_sea_water',
                'long_name': 'Silicate concentration in sea water',
                'units': 'umol L-1',
                'comment': 'Converted from original data assumed to be in mg/L to µmol/L',
                 'path': 'nutrients/silicate'},

    'amonio':{'v_name': 'ammonium',
               'nombre_fichero':f'{n_dataset}_NH4',
                'standard_name':'mole_concentration_of_ammonium_in_sea_water',
                'units': 'umol L-1',
                'long_name': 'Ammonium concentration in sea water',
                'comment': 'Converted from original data assumed to be in mg/L to µmol/L',
                 'path': 'nutrients/ammonium'},

    'NID':{'v_name': 'din',
            'nombre_fichero':f'{n_dataset}_DIN',
            'standard_name':'mole_concentration_of_dissolved_inorganic_nitrogen_in_sea_water',
            'long_name': 'Dissolved inorganic nitrogen concentration in sea water',
            'units': 'umol L-1',
            'comment': 'Converted from original data assumed to be in mg/L to µmol/L',
            'path': 'nutrients/din'},

    'relacion NP':{'v_name': 'NPrelation',
               'nombre_fichero':f'{n_dataset}_NPR',
                'standard_name':'',
                'long_name': 'Molar ratio of dissolved inorganic nitrogen to phosphate (N:P)',
                'units': 1,
                 'path': 'nutrients/relacionNP' },
}
data0 = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/bq/set9/resumen_belich_biogeoquimico.xlsx',  dtype={  
"Valor": "float64", }, parse_dates= ['Fecha'],  usecols= ['Variable', 'Estación', 'Fecha','Valor' ])
print(data0.head())
print(f'tiene una longitud de {len(data0)} filas')
print(data0.tail())
data0 = data0.replace(np.nan, -9999) #reemplazo los nan por -9999 
# %%
max_length_param = len("A")
estaciones = ['A', 'B', 'C', 'M']
estaciones_np = np.array(estaciones, dtype=f'S{max(len(s) for s in estaciones)}')
path = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/datasets_ncFormat/Biogeochemical/'
path_copia = 'C:/Users/Julia/Nextcloud/Datos_MM_Art_2025/Repository/Biogeochemical/'

# %%
def filtrar_df(data0, nombre_variable):
    data = data0.sort_values(by=["Fecha", "Estación"]).reset_index(drop=True)
    data  = data[data['Variable'] == nombre_variable]
    epoch = pd.Timestamp('1970-01-01')
    dias_desde_1970 = (data.Fecha.drop_duplicates() - epoch) / pd.Timedelta(days=1)
    print(data.head())
    print(f'tiene una longitud de {len(data)} filas')
    return data, dias_desde_1970
# %%
def crear_nc(data, dias_desde_1970, n_var):
    # nombre_variable_nc = nombre_variable.replace(":", "_")
    print(n_var)
    ncfile = Dataset(f"{path}{dict[n_var]['path']}/{n_dataset}/{dict[n_var]['nombre_fichero']}.nc", mode='w', format='NETCDF3_CLASSIC')

    ncfile.title =  dict[n_var]['nombre_fichero']
    ncfile.institution="Instituto Español de Oceanografía (IEO), Spain"
    ncfile.domain= 'Mar menor coastal lagoon, Spain'
    ncfile.dataset_id = n_dataset
    ncfile.project = 'BELICH'; ncfile.source = 'In situ data collection'; ncfile.Conventions = 'CF-1.8'

    ncfile.createDimension('time', len(dias_desde_1970))
    ncfile.createDimension('unit_char_len', max_length_param)
    ncfile.createDimension('station_name', len(estaciones))

    for dim in ncfile.dimensions.items():
        print(dim)

    time_var = ncfile.createVariable('time', np.float64, ('time',))
    time_var.units = "days since 1970-01-01 00:00:0"
    time_var.calendar = 'gregorian'
    time_var.standard_name = "time"
    time_var[:] = dias_desde_1970.values  # Se asigna directamente

    lat_var = ncfile.createVariable('latitude', np.float64,('station_name') )
    lat_var.units = 'degrees_north'
    lat_var.standard_name = "latitude"
    lat_var.grid_mapping = "crs"
    lat_var[:] = [ 37.791433,  37.70949,  37.667317, 37.710278, ]

    lon_var = ncfile.createVariable('longitude', np.float64,('station_name') )
    lon_var.units = 'degrees_east'
    lon_var.standard_name = "longitude"
    lon_var.grid_mapping = "crs"
    lon_var[:] = [ -0.78155 , -0.7851 , -0.755215, -0.831111]

    depth_var = ncfile.createVariable('depth',np.int8)
    depth_var.units = 'meters'
    depth_var.standard_name = 'depth'
    depth_var.positive = 'down'
    depth_var[:] = 4

    value_var = ncfile.createVariable(dict[n_var]['v_name'], np.float64, ('time', 'station_name'))
    value_var.units = dict[n_var]['units']
    value_var.standard_name = dict[n_var]['standard_name']
    value_var.long_name = dict[n_var]['long_name']
    value_var.missing_value = -9999
    value_var.grid_mapping = "crs"
    com = dict[n_var].get('comment')
    if com:
        value_var.comment = com

    # value_var.units = dict[n_var].get('units', '¿?¿')

    pivot = data.pivot_table(index='Fecha', columns='Estación', values='Valor')
    valor_array = pivot.to_numpy()
    value_var[:,:] = valor_array


    parameter_var = ncfile.createVariable('station', 'S1', ('station_name', 'unit_char_len'))
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

    ncfile.close()

# %% COMPROBACION
def comprobar_nc(n_var):
    # nombre_variable_nc = nombre_variable.replace(":", "_")
    dir = f"{path}{dict[n_var]['path']}/{n_dataset}/{dict[n_var]['nombre_fichero']}.nc"
    print(dir)
    dataset = Dataset(f"{path}{dict[n_var]['path']}/{n_dataset}/{dict[n_var]['nombre_fichero']}.nc", "r")
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
    unit = dataset.variables["station"][:]    #  
    value = dataset.variables[dict[nombre_variable]['v_name']][:] 

    print(tiempo); print('-----------------')
    print(unit); print('-----------------')
    print(value); print('-----------------')
    fechas = pd.to_datetime(tiempo, origin="1970-01-01", unit="D")
    fechas = np.array(fechas) 
    dataset.close()
    return fechas, value
# %% PLOT
def pintar_plot(fechas, value, nombre_variable, color = 'silver'):
    nombre_variable_nc = nombre_variable.replace(":", "_")
    fig, axs = plt.subplots(4, 1, figsize=(12, 10), sharex=True)

    for i, ax in enumerate(axs):
        ax.plot(fechas, value[:, i], label=nombre_variable, color=color, marker = 's')  
        ax.set_ylabel(nombre_variable)
        ax.set_title(f"Estación {estaciones[i]}")
        ax.legend()
        ax.grid(True)

    axs[-1].set_xlabel("Fecha")
    plt.tight_layout()
    plt.savefig(f'resumen_belich_bg_{nombre_variable_nc}.pdf', format='pdf')
    plt.show()

# %%

for nombre_variable in data0['Variable'].unique():
# nombre_variable= 'nitrato'
    # if nombre_variable == 'relacion N:P':
        data, dias_desde_1970 = filtrar_df(data0, nombre_variable)

        crear_nc(data, dias_desde_1970, nombre_variable)

        fechas, value= comprobar_nc(nombre_variable)

        # pintar_plot(fechas, value, nombre_variable, 'green')

        generar_txt(f"{path}{dict[nombre_variable]['path']}/{n_dataset}/{dict[nombre_variable]['nombre_fichero']}.nc", 
                   f"{path}{dict[nombre_variable]['path']}/{n_dataset}/{dict[nombre_variable]['nombre_fichero']}.txt")
        

        shutil.copy(f"{path}{dict[nombre_variable]['path']}/{n_dataset}/{dict[nombre_variable]['nombre_fichero']}.nc",
                    f"{path_copia}{dict[nombre_variable]['path']}/{n_dataset}/{dict[nombre_variable]['nombre_fichero']}.nc")

# %%