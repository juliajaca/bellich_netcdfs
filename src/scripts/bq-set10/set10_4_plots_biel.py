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
import matplotlib.pyplot as plt


%matplotlib widget
import matplotlib.pyplot as plt
# %%
data = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/bq/set10/sentinel_3_clorofila_BELA.xlsx',  dtype={  
"clorofila": "float64",}, parse_dates= ['fecha'], ) #rows=3
print(data.head())
print(f'tiene una longitud de {len(data)} filas')
# %%
fechas_unicas = np.sort(data["fecha"].unique())
fechas_unicas_ts = pd.to_datetime(fechas_unicas)
latitudes = np.sort(data["latitud"].unique())
longitudes = np.sort(data["longitud"].unique())

epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (fechas_unicas_ts - epoch) / pd.Timedelta(days=1)

# %%
fecha1 = data.loc[data['fecha'] == '2023-04-15 00:00:00']
print(f'la primera fecha tiene una longitud de {len(fecha1)} filas')

# %%
plt.figure(figsize=(8, 6))
sc1 = plt.scatter(fecha1.longitud, fecha1.latitud, cmap="viridis", marker="o", edgecolors="k", label="7 enero 2023", alpha= 0.5, s= 10)
# sc2 = plt.scatter(fecha2.longitud, fecha2.latitud, c=fecha2.clorofila, cmap="coolwarm", marker="s", edgecolors="k", label="11 enero2023 ", alpha = 0.5, s=10)
# sc3 = plt.scatter(fecha3.longitud, fecha3.latitud, c=fecha3.clorofila, cmap="jet", marker="x", edgecolors="k", label="10 enero 2025", alpha = 0.5, s=10)
# plt.colorbar(sc1, label="Clorofila")
plt.xlabel("Longitud")
plt.ylabel("Latitud")
plt.title("Clorofila para el 7 y el 11 de enero")
# plt.legend()
plt.grid()
plt.show()

# %%
plt.figure(figsize=(12, 8))
sc1 = plt.scatter(fecha1[0:100].longitud, fecha1[0:100].latitud, cmap="viridis", marker="o", edgecolors="k", label="7 enero 2023", alpha= 0.5, s= 10)

for i in range(100):
    x = fecha1.iloc[i].longitud
    y = fecha1.iloc[i].latitud
    plt.text(x + 0.0001, y + 0.0001, str(fecha1.index[i]), fontsize=16)
plt.xlabel("Longitud"); plt.ylabel("Latitud"); plt.title("Clorofila para el 7 de enero");plt.grid();plt.show()
# %%
# Índices que quieres destacar
destacados = [0, 1, 2, 3,4,5]

# Pintar puntos normales (en negro)
for i in range(len(fecha1)):
    if i not in destacados:
        plt.scatter(fecha1.iloc[i].longitud, fecha1.iloc[i].latitud, color='black', s=10, edgecolors='k', alpha=0.5)
    else:
        plt.scatter(fecha1.iloc[i].longitud, fecha1.iloc[i].latitud, color='red', s=80, edgecolors='k', alpha=0.5)
        plt.text(fecha1.iloc[i].longitud + 0.001, fecha1.iloc[i].latitud + 0.001, str(fecha1.index[i]), fontsize=8)

# Pintar puntos destacados (en rojo y más grandes)
# for i in destacados:
#     plt.scatter(longitudes[i], latitudes[i], color='red', s=50, edgecolors='k', alpha=0.8)

# Etiquetas
# for i in (destacados):
#     plt.text(longitudes[i] + 0.001, latitudes[i] + 0.001, str(fecha1.index[i]), fontsize=8)

plt.xlabel("Longitud")
plt.ylabel("Latitud")
plt.title("Clorofila para el 7 de enero")
plt.grid()
plt.show()
# %% ver lat longs
for i in range(100):
    print('punto',i, ': ',fecha1.iloc[i].longitud, ', ',fecha1.iloc[i].latitud)
# %%
# ---------------------------------------------------------
distancias_lat = pd.Series( [(fecha1.iloc[i+1].latitud - fecha1.iloc[i].latitud) for i in range(len(fecha1) - 1)])
distancias_lon = pd.Series([(fecha1.iloc[i+1].longitud - fecha1.iloc[i].longitud) for i in range(len(fecha1) - 1)])
mediana_lat = distancias_lat.median() * 1.5
mediana_lon = distancias_lon.median() * 1.5
saltos = (np.abs(distancias_lat) > np.abs(mediana_lat)) | (np.abs(distancias_lon) > np.abs(mediana_lon))

# %%
plt.figure(figsize=(10, 5))
plt.axhline( distancias_lon.median(), color='blue', linestyle='--', label=f'Mediana dx: {distancias_lon.median():.2e}')
plt.axhline(distancias_lat.median(), color='orange', linestyle='--', label=f'Mediana dy: {distancias_lat.median():.2e}')
plt.plot(distancias_lon, label='Δ Longitud (dx)', marker='o')
plt.plot(distancias_lat, label='Δ Latitud (dy)', marker='x')
plt.title('Diferencias entre puntos consecutivos')
plt.xlabel('Índice')
plt.ylabel('Diferencia')
plt.legend()
plt.grid(True)
plt.show()

# %%
lat_nueva = []
lon_nueva = []
chl_nueva = []

for i in range(len(fecha1) - 1):
    # Añadimos el punto actual
    lat_nueva.append(fecha1.iloc[i]['latitud'])
    lon_nueva.append(fecha1.iloc[i]['longitud'])
    chl_nueva.append(fecha1.iloc[i]['clorofila'])

    # Si hay salto, añadimos un NaN
    if saltos[i]:
        lat_nueva.append(fecha1.iloc[i]['latitud']+ distancias_lat.median())
        lon_nueva.append(fecha1.iloc[i]['longitud']+ distancias_lon.median())
        chl_nueva.append(np.nan)
        lat_nueva.append(fecha1.iloc[i+1]['latitud']- distancias_lat.median())
        lon_nueva.append(fecha1.iloc[i+1]['longitud']- distancias_lon.median())
        chl_nueva.append(np.nan)

# Añadir el último punto
lat_nueva.append(fecha1.iloc[len(fecha1)-1]['latitud'])
lon_nueva.append(fecha1.iloc[len(fecha1)-1]['longitud'])
chl_nueva.append(fecha1.iloc[len(fecha1)-1]['clorofila'])

# Construir nuevo DataFrame con los NaNs incluidos
df_completo = pd.DataFrame({'latitud': lat_nueva,'longitud': lon_nueva,'clorofila': chl_nueva})
print(df_completo.head(15))  # vista previa

# %%
# Filtrar puntos donde 'clorofila' es NaN
mask_nan = df_completo['clorofila'].isna()

plt.figure(figsize=(12, 8))

plt.scatter(df_completo.loc[mask_nan, 'longitud'], df_completo.loc[mask_nan, 'latitud'], c='red', label='Clorofila NaN', s=40, edgecolors='black', alpha=0.5)
plt.scatter(df_completo.loc[~mask_nan, 'longitud'], df_completo.loc[~mask_nan, 'latitud'], c='blue', label='Valores válidos', s=20)
# for i in range((2000)):
#     x = df_completo.iloc[i].longitud
#     y = df_completo.iloc[i].latitud
#     plt.text(x + 0.0001, y + 0.0001, str(df_completo.index[i]), fontsize=16)
plt.xlabel('Longitud')
plt.ylabel('Latitud')
plt.title('Distribución de Clorofila metodo biel')

plt.legend()
plt.grid(True)
plt.show()
# %%
lc = pd.read_csv('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/bq/set10/marmenor_inside.dat', delim_whitespace=True, header=None, names=['x', 'y', 'z'])
print(lc.head())
lc2 = pd.read_csv('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/bq/set10/marmenor_coastline.dat', delim_whitespace=True, header=None, names=['x', 'y', 'z'])
print(lc2.head())
bat =pd.concat([lc2,lc])
print(bat.head())
# %%
plt.figure(figsize=(6, 6))
plt.scatter(lc['x'], lc['y'],  label='Puntos', s= 10)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Plot de posiciones')
plt.grid(True)
plt.legend()
plt.axis('equal')  # para que los ejes tengan la misma escala
# for i, (xi, yi) in enumerate(zip(lc['x'], lc['y'])):
#     plt.text(xi + 0.005, yi + 0.005, str(i), fontsize=8, color='black')
plt.show()
# %%
#  CREAR LA MALLA
# Número de puntos en cada dirección
n_lat = 50
n_lon = 50
n_lat = int(np.ceil((latitudes.max() - latitudes.min()) / abs(distancias_lat.median())))
# n_lon = int(np.ceil((longitudes.max() - longitudes.min()) / distancias_lon.median()))

n_lon = int(np.ceil((longitudes.max() - longitudes.min()) / abs(distancias_lat.median())))
# n_lat = int(np.ceil((latitudes.max() - latitudes.min()) / abs(distancias_lon.median())))

# Crear la malla
latitudes = np.linspace(latitudes.min(), latitudes.max(), n_lat)
longitudes = np.linspace(longitudes.min(), longitudes.max(), n_lon)
lon_grid, lat_grid = np.meshgrid(longitudes, latitudes)

# %%
# Mostrar la malla en un gráfico
plt.figure(figsize=(8, 6))
plt.scatter(lon_grid, lat_grid, s=10, color='blue')
plt.scatter(bat['x'], bat['y'],  s=10, color='red', label='Batimetría')
plt.title("Malla 2D de coordenadas")
plt.xlabel("Longitud")
plt.ylabel("Latitud")
plt.grid(True)
plt.axis('equal')
plt.show()


# %%
# Inteprolañar
bat = bat.dropna(subset=['x', 'y', 'z']) #qiuto los nans para interpolar
puntos = np.column_stack((bat['x'], bat['y']))
valores = bat['z'].values

# Interpolación sobre la malla
batimetria_interpolada = griddata(puntos, valores, (lon_grid, lat_grid), method='linear')

# %%
plt.figure(figsize=(10, 8))
plt.contourf(lon_grid, lat_grid, batimetria_interpolada, cmap='viridis', levels=100)
plt.colorbar(label='Profundidad interpolada (m)')
plt.title("Mapa Interpolado de Batimetría")
plt.xlabel("Longitud")
plt.ylabel("Latitud")
plt.grid(True)
plt.axis('equal')
plt.show()
# %% CREAR MASCARA
mascara = (batimetria_interpolada <= 0) & (batimetria_interpolada >= -1) | ((np.isnan(batimetria_interpolada)))

# %% INTERPOlAR FECHA 1
punts_fecha1 = np.column_stack((df_completo.longitud, df_completo.latitud))
valores_fecha1 = df_completo.clorofila.values
# %%
fecha1_interpolada = griddata(punts_fecha1, valores_fecha1, (lon_grid, lat_grid), method='linear')
print(f' hay {(~np.isnan(fecha1_interpolada)).sum()} valores no nan en la interpolacion de la fecha')
print(f' hay {(~np.isnan(fecha1.clorofila)).sum()} valores no nan en la  fecha real')
fecha1_interpolada = griddata(punts_fecha1, valores_fecha1, (lon_grid, lat_grid), method='cubic')
print(f' hay {(~np.isnan(fecha1_interpolada)).sum()} valores no nan en la interpolacion de la fecha')
print(f' hay {(~np.isnan(fecha1.clorofila)).sum()} valores no nan en la  fecha real')
fecha1_interpolada = griddata(punts_fecha1, valores_fecha1, (lon_grid, lat_grid), method='nearest')
print(f' hay {(~np.isnan(fecha1_interpolada)).sum()} valores no nan en la interpolacion de la fecha')
print(f' hay {(~np.isnan(fecha1.clorofila)).sum()} valores no nan en la  fecha real')
# Interpolación lineal
fecha1_interpolada = griddata(punts_fecha1, valores_fecha1, (lon_grid, lat_grid), method='linear')
# Dónde hay NaN
mask_nan = np.isnan(fecha1_interpolada)
# Interpolar solo esos puntos usando nearest
fecha1_interpolada[mask_nan] = griddata(punts_fecha1, valores_fecha1, (lon_grid[mask_nan], lat_grid[mask_nan]), method='nearest')
print(f' hay {(~np.isnan(fecha1_interpolada)).sum()} valores no nan en la interpolacion de la fecha')
print(f' hay {(~np.isnan(fecha1.clorofila)).sum()} valores no nan en la  fecha real')
# %%
plt.figure(figsize=(10, 8))
plt.scatter(lon_grid, lat_grid, c= fecha1_interpolada, cmap='viridis',s=100)
plt.colorbar(label='Clorofila interpolada ')
plt.title("Mapa Interpolado de Batimetría")
plt.xlabel("Longitud")
plt.ylabel("Latitud")
plt.grid(True)
plt.axis('equal')
plt.show()

# %%
# APLICAR mascara
# Paso 1: Obtener las coordenadas de longitud y latitud donde la máscara es True
lon_mascara = lon_grid[~mascara]
lat_mascara = lat_grid[~mascara]

# Paso 2: Obtener los valores de clorofila correspondientes a esos puntos
clorofila_mascarada = fecha1_interpolada[~mascara]

# Paso 3: Crear el DataFrame con latitud, longitud y clorofila
df_clorofila_mascarada = pd.DataFrame({
    'longitud': lon_mascara,
    'latitud': lat_mascara,
    'clorofila': clorofila_mascarada
})

# Mostrar las primeras filas del DataFrame
print(df_clorofila_mascarada.head())
# %%
plt.figure(figsize=(10, 8))
plt.scatter(df_clorofila_mascarada.longitud, df_clorofila_mascarada.latitud, c=df_clorofila_mascarada.clorofila)
plt.colorbar(label='Concentración de Clorofila (mg/m³)')
plt.title('Clorofila Interpolada en la Zona de Batimetría entre 0 y -1 m')
plt.xlabel('Longitud')
plt.ylabel('Latitud')
plt.legend()
plt.grid(True)
plt.axis('equal')
plt.show()

# %%  PLOT BIEL SUPERPONER MASCARA + SATELITE
plt.figure(figsize=(8, 6))
plt.scatter(data.longitud, data.latitud, color='black', s=10, edgecolors='k', alpha=0.5,label='datos reales')
# plt.scatter(fecha1.longitud, fecha1.latitud, color='black', s=10, edgecolors='k', alpha=0.5,label='datos reales fecha1')
plt.scatter(df_clorofila_mascarada.longitud, df_clorofila_mascarada.latitud, color = 'red', alpha = 0.5, s = 10, label='máscara')
plt.title("SUperponer mascara + satelite")
plt.xlabel("Longitud")
plt.ylabel("Latitud")
plt.grid(True)
plt.legend() 
plt.axis('equal')
plt.show()

# %% plot biel 2 paneles --> interpolacion + satelite
fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True)
sc1=axes[0].scatter(fecha1.longitud, fecha1.latitud, c= fecha1.clorofila, label='datos reales')
axes[0].set_title("Datos reales fecha 1")
axes[0].set_xlabel("Longitud")
axes[0].set_ylabel("Latitud")
axes[0].grid(True)
axes[0].legend()
fig.colorbar(sc1, ax=axes[0], label='Clorofila real')
sc2=axes[1].scatter(lon_grid, lat_grid, c=fecha1_interpolada,  label='fecha 1 interpolada')
axes[1].set_title("fehca 1 interpolada ")
axes[1].set_xlabel("Longitud")
axes[1].grid(True)
axes[1].legend()
fig.colorbar(sc2, ax=axes[1], label='Clorofila Interpolada')
plt.tight_layout()
plt.show()

# %% MAPA INTERPOLADO DE LA BATIMETRIA
plt.figure(figsize=(10, 8))
plt.contourf(lon_grid, lat_grid, batimetria_interpolada, cmap='viridis', levels=100)
plt.colorbar(label='Clorofila interpolada ')
plt.title("Mapa Interpolado de Batimetría")
plt.xlabel("Longitud")
plt.ylabel("Latitud")
plt.grid(True)

plt.show()
# %%
plt.figure(figsize=(10, 8))
plt.contourf(lon_grid, lat_grid, mascara, cmap='viridis', levels=100)
plt.colorbar(label='Clorofila interpolada ')
plt.title("Mapa Interpolado de Batimetría")
plt.xlabel("Longitud")
plt.ylabel("Latitud")
plt.grid(True)

plt.show()

# %%
