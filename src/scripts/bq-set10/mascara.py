# %%
from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import geopandas as gpd
import alphashape
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from geopy.distance import geodesic
from matplotlib.path import Path


# %%
# %%
def get_mascara(lon_grid, lat_grid):
    lc = pd.read_csv('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/bq/set10/marmenor_inside.dat', delim_whitespace=True, header=None, names=['x', 'y', 'z'])
    print(lc.head())
    lc2 = pd.read_csv('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/bq/set10/marmenor_coastline.dat', delim_whitespace=True, header=None, names=['x', 'y', 'z'])
    print(lc2.head())
    bat = pd.concat([lc2,lc])
    bat = bat.dropna(subset=['x', 'y', 'z'])
    puntos = np.column_stack((bat['x'], bat['y']))
    valores = bat['z'].values
    # Interpolación sobre la malla
    batimetria_interpolada = griddata(puntos, valores, (lon_grid, lat_grid), method='linear')
    mascara = (batimetria_interpolada <= 0) & (batimetria_interpolada >= -1) | ((np.isnan(batimetria_interpolada)))
    return mascara

def get_df_clorofila_mascarada(lon_grid, lat_grid, mascara, clorofila_interpolada, fecha):
    # Obtener la longitud,latitud y chl donde la máscara es True
    lon_mascara = lon_grid[~mascara]
    lat_mascara = lat_grid[~mascara]
    clorofila_mascara = clorofila_interpolada[~mascara]
    # Crear el DataFrame con latitud, longitud y clorofila
    df_clorofila_mascarada = pd.DataFrame({
        'fecha': fecha,
        'longitud': lon_mascara,
        'latitud': lat_mascara,
        'clorofila': clorofila_mascara
    })
    return df_clorofila_mascarada