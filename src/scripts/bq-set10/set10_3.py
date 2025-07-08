# %%
from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from geopy.distance import geodesic
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
# %%
# Calcular la diferencia en días
dias_desde_1970 = (fechas_unicas_ts - epoch) / pd.Timedelta(days=1)

# %%
fecha1 = data.loc[data['fecha'] == '2023-01-07 00:00:00']
print(f'la primera fecha tiene una longitud de {len(fecha1)} filas')
fecha2 = data.loc[data['fecha'] == '2023-01-11T00:00:00.000000000'] 
print(f'la segunda fecha tiene una longitud de {len(fecha2)} filas')
fecha3 = data.loc[data['fecha'] == '2025-01-10T00:00:00.000000000'] 
print(f'la segunda fecha tiene una longitud de {len(fecha3)} filas')
fecha4 = data.loc[data['fecha'] == '2023-07-13T00:00:00.000000000'] 
print(f'la segunda fecha tiene una longitud de {len(fecha4)} filas')


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

plt.xlabel("Longitud")
plt.ylabel("Latitud")
plt.title("Clorofila para el 7 de enero")

plt.grid()

plt.show()
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


# %% VECTORIZADO
# # Paso 1: Obtener las coordenadas como lista de tuplas (lat, lon)
# coords = list(zip(fecha1['latitud'], fecha1['longitud']))

# # Paso 2: Calcular las distancias entre puntos consecutivos (vectorizado usando una lista de comprensión)
# distancias = np.array([geodesic(coords[i], coords[i+1]).meters for i in range(len(coords) - 1)])

# # Paso 3: Calcular el umbral de saltos inusuales
# umbral = np.median(distancias) * 1.5

# # Paso 4: Detectar los saltos inusuales (crear un array booleano)
# saltos = distancias > umbral

# # Paso 5: Crear las listas de latitud, longitud y clorofila
# # Asegúrate de que las dimensiones de fecha1 sean compatibles con el número de puntos a procesar
# lat_nueva = np.repeat(fecha1['latitud'].values, 2)[:-1]  # Repetimos las latitudes para cada par de puntos
# lon_nueva = np.repeat(fecha1['longitud'].values, 2)[:-1]  # Lo mismo para longitudes
# chl_nueva = np.repeat(fecha1['clorofila'].values, 2)[:-1]  # Y para clorofila

# # Paso 6: Detectar dónde hay saltos y reemplazar con NaN
# lat_nueva[saltos] = np.nan
# lon_nueva[saltos] = np.nan
# chl_nueva[saltos] = np.nan

# # Paso 7: Añadir el último punto (que no se repitió en el paso 5)
# lat_nueva = np.append(lat_nueva, fecha1.iloc[-1]['latitud'])
# lon_nueva = np.append(lon_nueva, fecha1.iloc[-1]['longitud'])
# chl_nueva = np.append(chl_nueva, fecha1.iloc[-1]['clorofila'])

# # Paso 8: Construir el nuevo DataFrame con los NaNs incluidos
# df_completo = pd.DataFrame({
#     'latitud': lat_nueva,
#     'longitud': lon_nueva,
#     'clorofila': chl_nueva
# })

# print(df_completo.head(15))  # Vista previa
# ---------------------------------------------------------
# Paso 1: calcular distancia entre puntos consecutivos
coords = list(zip(fecha1['latitud'], fecha1['longitud']))
distancias = [geodesic(coords[i], coords[i+1]).meters for i in range(len(coords) - 1)]
distancias = pd.Series(distancias)

# Paso 2: detectar saltos inusuales (umbral = 1.5 * mediana)
umbral = distancias.median() * 1.5
saltos = distancias > umbral

# Paso 3: construir nuevo DataFrame insertando NaNs donde hay huecos
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
        lat_nueva.append(np.nan)
        lon_nueva.append(np.nan)
        chl_nueva.append(np.nan)

# Añadir el último punto
lat_nueva.append(fecha1.iloc[len(fecha1)-1]['latitud'])
lon_nueva.append(fecha1.iloc[len(fecha1)-1]['longitud'])
chl_nueva.append(fecha1.iloc[len(fecha1)-1]['clorofila'])

# Construir nuevo DataFrame con los NaNs incluidos
df_completo = pd.DataFrame({
    'latitud': lat_nueva,
    'longitud': lon_nueva,
    'clorofila': chl_nueva
})

print(df_completo.head(15))  # vista previa


# %%
lats_def = []
lons_def= []
chls_def= []
grupo_lat = []
grupo_lon =[]
grupo_chl = []
for fila in range(len(df_completo[0:])):
    # print(df_completo.loc[fila])
    if pd.isna(df_completo.loc[fila].longitud):
        print('es nan!')
        print(grupo_lat)
        print(grupo_lon)
        lat_direccion = grupo_lat[-1] - grupo_lat[0]
        lon_direccion = grupo_lon[-1] - grupo_lon[0]
        lat_nueva_antes = grupo_lat[0] - lat_direccion * 0.1
        lat_nueva_despues= grupo_lat[-1] + lat_direccion* 0.1
        lon_nueva_antes= grupo_lon[0]  - lon_direccion* 0.1
        lon_nueva_despues = grupo_lon[-1]  + lon_direccion* 0.1
        print(f'la lat direccion es {lat_direccion} y la lon {lon_direccion} ')
        print(f'posicion nueva antes {lat_nueva_antes} y {lon_nueva_antes}')
        print(f'posicion nueva despues {lat_nueva_despues} y {lon_nueva_despues}')
        lats_def.extend([lat_nueva_antes])
        lats_def.extend(grupo_lat)
        lats_def.extend([lat_nueva_despues])

        lons_def.extend([lon_nueva_antes])
        lons_def.extend(grupo_lon)
        lons_def.extend([lon_nueva_despues])

        chls_def.extend([np.nan])
        chls_def.extend(grupo_chl)
        chls_def.extend([np.nan])


        #vacio los grupos
        grupo_lat = []
        grupo_lon =[]
        grupo_chl = []
    else:
        # print(df_completo.loc[fila])
        grupo_lat.append(df_completo.loc[fila].latitud)
        grupo_lon.append(df_completo.loc[fila].longitud)
        grupo_chl.append(df_completo.loc[fila].clorofila)

    print('---------------------')

df_final = pd.DataFrame({
    'latitud': lats_def,
    'longitud': lons_def,
    'clorofila': chls_def
})

# %% PLOT
# Filtrar puntos donde 'clorofila' es NaN
mask_nan = df_final['clorofila'].isna()

# Crear el gráfico de dispersión
plt.figure(figsize=(12, 8))

# Puntos con clorofila NaN (rojos y más grandes)
plt.scatter(df_final.loc[mask_nan, 'longitud'], df_final.loc[mask_nan, 'latitud'], c='red', label='Clorofila NaN', s=40, edgecolors='black', alpha=0.5)
# Puntos con valores válidos de clorofila
plt.scatter(df_final.loc[~mask_nan, 'longitud'], df_final.loc[~mask_nan, 'latitud'], c='blue', label='Valores válidos', s=20)


# Etiquetas y título
plt.xlabel('Longitud')
plt.ylabel('Latitud')
plt.title('Distribución de Clorofila')

# Mostrar leyenda
plt.legend()

# Mostrar el gráfico
plt.grid(True)
plt.show()
# %%
