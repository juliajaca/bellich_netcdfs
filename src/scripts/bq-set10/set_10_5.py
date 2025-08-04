# %%
from netCDF4 import Dataset, date2num,num2date, stringtochar
import numpy as np
import geopandas as gpd
import matplotlib.gridspec as gridspec
import alphashape
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.spatial import cKDTree
from geopy.distance import geodesic
from matplotlib.path import Path
import matplotlib.pyplot as plt
from mascara import get_df_clorofila_mascarada, get_mascara
# %%%
def pintar_destacados(fecha1):
    destacados = [0, 1, 2, 3,4,5]

    # Pintar puntos normales (en negro)
    for i in range(len(fecha1)):
        if i not in destacados:
            plt.scatter(fecha1.iloc[i].longitud, fecha1.iloc[i].latitud, color='black', s=10, edgecolors='k', alpha=0.5)
        else:
            plt.scatter(fecha1.iloc[i].longitud, fecha1.iloc[i].latitud, color='red', s=80, edgecolors='k', alpha=0.5)
            plt.text(fecha1.iloc[i].longitud + 0.001, fecha1.iloc[i].latitud + 0.001, str(fecha1.index[i]), fontsize=8)

def hacer_plot_nan(df_completo):
    mask_nan = df_completo['clorofila'].isna()
    plt.figure(figsize=(12, 8))
    plt.scatter(df_completo.loc[mask_nan, 'longitud'], df_completo.loc[mask_nan, 'latitud'], c='red', label='Clorofila NaN', s=70, edgecolors='black', alpha=0.5)
    plt.scatter(df_completo.loc[~mask_nan, 'longitud'], df_completo.loc[~mask_nan, 'latitud'], c='blue', label='Valores válidos', s=20)
    plt.xlabel('Longitud')
    plt.ylabel('Latitud')
    plt.title('Distribución de Clorofila metodo biel')
    plt.legend()
    plt.grid(True)
    plt.show()
# %%
def hacer_plot(df_clorofila_mascarada):
    # Encuentra el rango común
    vmin = min(fecha1.clorofila.min(), df_clorofila_mascarada.clorofila.min())
    vmax = max(fecha1.clorofila.max(), df_clorofila_mascarada.clorofila.max())

    # Ahora pasa vmin y vmax a cada scatter:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True)

    sc1 = axes[0].scatter(fecha1.longitud, fecha1.latitud, c=fecha1.clorofila, vmin=vmin, vmax=vmax, )
    axes[0].contourf(lon_grid, lat_grid, mascara, cmap='viridis', levels=2, alpha = 0.1 )
    # axes[0].scatter(lon_grid, lat_grid, s= 0.01 )
    axes[0].set_title(f"Datos reales fecha {fecha}")
    axes[0].set_xlabel("Longitud")
    axes[0].set_ylabel("Latitud")
    axes[0].grid(True)
    axes[0].text(0.02, 0.02, f'nº de puntos: {len(fecha1)}', transform=axes[0].transAxes, fontsize=10, ha='left', va='bottom')
    fig.colorbar(sc1, ax=axes[0], label='Clorofila real')

    sc2 = axes[1].scatter(df_clorofila_mascarada.longitud, df_clorofila_mascarada.latitud, c=df_clorofila_mascarada.clorofila, vmin=vmin, vmax=vmax,  )
    axes[1].set_title(f"Datos interpolados fecha {fecha}")
    axes[1].contourf(lon_grid, lat_grid, mascara, cmap='viridis', levels=2, alpha = 0.1 )
    axes[1].set_xlabel("Longitud")
    axes[1].grid(True)
    axes[1].text(0.02, 0.02, f'nº de puntos: {(~np.isnan(df_clorofila_mascarada.clorofila)).sum()}', transform=axes[1].transAxes, fontsize=10, ha='left', va='bottom')
    fig.colorbar(sc2, ax=axes[1], label='Clorofila Interpolada y mascarada')
    plt.tight_layout()

    # plt.figure(figsize=(8, 6))
    # plt.scatter(data.longitud, data.latitud, color='black', s=10, edgecolors='k', alpha=0.5,label='datos reales')
    # # plt.scatter(fecha1.longitud, fecha1.latitud, color='black', s=10, edgecolors='k', alpha=0.5,label='datos reales fecha1')
    # plt.scatter(df_clorofila_mascarada.longitud, df_clorofila_mascarada.latitud, color = 'red', alpha = 0.5, s = 10, label='máscara')
    # plt.title("SUperponer mascara + satelite")
    # plt.xlabel("Longitud")
    # plt.ylabel("Latitud")
    # plt.grid(True)
    # plt.legend() 
    # plt.axis('equal')
    # plt.show()
    timestamp = fecha.strftime('%Y%m%d')
    plt.savefig(f'2_sentinel_3_clorofila_bela_{timestamp}.pdf', format='pdf')
    plt.show()

def hacer_5_plots(df1, df2,df3,df4,df5):
    vmin = min(df1.clorofila.min(), df2.clorofila.min(), df3.clorofila.min(),df4.clorofila.min(), df5.clorofila.min())
    vmax = max(df1.clorofila.max(), df2.clorofila.max(),df3.clorofila.max(), df4.clorofila.max(), df5.clorofila.max())
    # Crear figura y layout: 3 filas x 2 columnas
    fig = plt.figure(figsize=(8.5, 11))  # Tamaño vertical A4
    gs = gridspec.GridSpec(3, 2)

    # Primer gráfico: fila 0, columna 0
    ax1 = fig.add_subplot(gs[0, 0])
    sc1 = ax1.scatter(df1.longitud, df1.latitud, c=df1.clorofila, vmin=vmin, vmax=vmax, s=5)
    ax1.contourf(lon_grid, lat_grid, mascara, cmap='viridis', levels=2, alpha=0.1)
    ax1.set_title(f"Datos reales fecha {fecha}")
    ax1.set_xlabel("Longitud")
    ax1.set_ylabel("Latitud")
    ax1.grid(True)
    ax1.text(0.02, 0.02, f'nº de puntos: {len(df1)}', transform=ax1.transAxes, fontsize=10, ha='left', va='bottom')
    fig.colorbar(sc1, ax=ax1, label='Clorofila real')

    # Segundo gráfico: fila 0, columna 1
    ax2 = fig.add_subplot(gs[0, 1])
    sc2 = ax2.scatter(df2.longitud, df2.latitud, 
                    c=df2.clorofila, vmin=vmin, vmax=vmax, s=5)
    ax2.contourf(lon_grid, lat_grid, mascara, cmap='viridis', levels=2, alpha=0.1)
    ax2.set_title(f"Interpolado {fecha} linear")
    ax2.set_xlabel("Longitud")
    ax2.set_ylabel("Latitud")
    ax2.grid(True)
    ax2.text(0.02, 0.02, f'nº de puntos: {(~np.isnan(df2.clorofila)).sum()}',transform=ax2.transAxes, fontsize=10, ha='left', va='bottom')
    fig.colorbar(sc2, ax=ax2, label='Interpolada linear')

    # Los otros tres slots están listos para más plots
    ax3 = fig.add_subplot(gs[1, 0])
    sc3 = ax3.scatter(df3.longitud, df3.latitud, 
                    c=df3.clorofila, vmin=vmin, vmax=vmax, s=5)
    ax3.contourf(lon_grid, lat_grid, mascara, cmap='viridis', levels=2, alpha=0.1)
    ax3.set_title(f"Interpolado {fecha} nearest")
    ax3.set_xlabel("Longitud")
    ax3.set_ylabel("Latitud")
    ax3.grid(True)
    ax3.text(0.02, 0.02, f'nº de puntos: {(~np.isnan(df3.clorofila)).sum()}',transform=ax3.transAxes, fontsize=10, ha='left', va='bottom')
    fig.colorbar(sc3, ax=ax3, label='Interpolada nearest')

    ax4 = fig.add_subplot(gs[1, 1])
    sc4  = ax4.scatter(df4.longitud, df4.latitud, 
                    c=df4.clorofila, vmin=vmin, vmax=vmax, s=5)
    ax4.contourf(lon_grid, lat_grid, mascara, cmap='viridis', levels=2, alpha=0.1)
    ax4.set_title(f"Interpolado {fecha} linear nearest")
    ax4.set_xlabel("Longitud")
    ax4.set_ylabel("Latitud")
    ax4.grid(True)
    ax4.text(0.02, 0.02, f'nº de puntos: {(~np.isnan(df4.clorofila)).sum()}',transform=ax4.transAxes, fontsize=10, ha='left', va='bottom')
    fig.colorbar(sc4, ax=ax4, label='Interpolada linear nearest')

    ax5 = fig.add_subplot(gs[2, 0])
    sc5  = ax5.scatter(df5.longitud, df5.latitud, 
                    c=df5.clorofila, vmin=vmin, vmax=vmax, s=5)
    ax5.contourf(lon_grid, lat_grid, mascara, cmap='viridis', levels=2, alpha=0.1)
    ax5.set_title(f"Interpolado {fecha} nearest cerca")
    ax5.set_xlabel("Longitud")
    ax5.set_ylabel("Latitud")
    ax5.grid(True)
    ax5.text(0.02, 0.02, f'nº de puntos: {(~np.isnan(df5.clorofila)).sum()}',transform=ax5.transAxes, fontsize=10, ha='left', va='bottom')
    fig.colorbar(sc5, ax=ax5, label='Interpolada  nearest cerca')

    # Si no quieres usar gs[2,1], podrías dejarlo en blanco o eliminarlo con:
    fig.delaxes(fig.add_subplot(gs[2, 1]))

    plt.tight_layout()
    timestamp = fecha.strftime('%Y%m%d')
    # plt.savefig(f'sentinel_3_clorofila_bela_{timestamp}_distancias_lat.medianX1.5.pdf', format='pdf')
    plt.show()
# %%

# lc = pd.read_csv('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/bq/set10/marmenor_inside.dat', delim_whitespace=True, header=None, names=['x', 'y', 'z'])
# print(lc.head())
# lc2 = pd.read_csv('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/bq/set10/marmenor_coastline.dat', delim_whitespace=True, header=None, names=['x', 'y', 'z'])
# print(lc2.head())
# bat =pd.concat([lc2,lc])
# bat = lc2[
#     (lc2['x'] <= -0.72)  & 
#     (lc2['y'] >= 37.63) & (lc2['y'] <= 37.82)
# ]

data = pd.read_excel('C:/Users/Julia/Documents/VSCODE_BELLICH/src/datos/bq/set10/sentinel_3_clorofila_BELA.xlsx',  dtype={  
"clorofila": "float64",}, parse_dates= ['fecha'], ) #rows=3
print(data.head())
print(f'tiene una longitud de {len(data)} filas')
# %%
fechas_unicas = np.sort(data["fecha"].unique())
fechas_unicas_ts = pd.to_datetime(fechas_unicas)
epoch = pd.Timestamp('1970-01-01')
dias_desde_1970 = (fechas_unicas_ts - epoch) / pd.Timedelta(days=1)

# %% MALLA
fecha1 = data.loc[data['fecha'] == '2023-01-07 00:00:00']

latitudes = np.sort(data["latitud"].unique())
longitudes = np.sort(data["longitud"].unique())

distancias_lat = pd.Series( [(fecha1.iloc[i+1].latitud - fecha1.iloc[i].latitud) for i in range(len(fecha1) - 1)])
distancias_lon = pd.Series([(fecha1.iloc[i+1].longitud - fecha1.iloc[i].longitud) for i in range(len(fecha1) - 1)])
mediana_lat = distancias_lat.median() * 1.5
mediana_lon = distancias_lon.median() * 1.5

n_lat = int(np.ceil((latitudes.max() - latitudes.min()) / abs(distancias_lat.median())))
n_lon = int(np.ceil((longitudes.max() - longitudes.min()) / abs(distancias_lat.median())))
# n_lat = 55
# n_lon = 41
# creamos malla
latitudes = np.linspace(latitudes.min(), latitudes.max(), n_lat)
longitudes = np.linspace(longitudes.min(), longitudes.max(), n_lon)
lon_grid, lat_grid = np.meshgrid(longitudes, latitudes)
mascara = get_mascara(lon_grid, lat_grid)

# %%
df_final = pd.DataFrame()
fechas_malas = ['20230228']
# fechas_malas =  ['20230228',
# '20230304','20230323','20230327','20230331','20230516','20230713','20230508','20230928','20230423',
#                 '20231002', '20231006','20231029','20231102','20231125','20231129','20231218','20231222','20240702','20240721',
#                 '20240817', '20241018','20241110','20241129','2024127',]
for fecha in data.fecha.unique()[0:]:
    # print(fecha.strftime('%Y%m%d'))
    # if fecha.strftime('%Y%m%d') in fechas_malas:
        fecha1 = data.loc[data['fecha'] == fecha].reset_index()
        # print(fecha1.head(10))
        fecha1 = fecha1.drop_duplicates(subset=["latitud", "longitud"])
        print('----')
        # print(fecha1.head(10))
        print(f'la fecha {fecha} tiene una longitud de {len(fecha1)} filas')

        # pintar_destacados(fecha1)
        distancias_lat = pd.Series( [(fecha1.iloc[i+1].latitud - fecha1.iloc[i].latitud) for i in range(len(fecha1) - 1)])
        distancias_lon = pd.Series([(fecha1.iloc[i+1].longitud - fecha1.iloc[i].longitud) for i in range(len(fecha1) - 1)])

        mediana_lat = distancias_lat.median() * 1.5
        mediana_lon = distancias_lon.median() * 1.5
        saltos = (np.abs(distancias_lat) > np.abs(mediana_lat)) | (np.abs(distancias_lon) > np.abs(mediana_lon))
        lat_nueva = []
        lon_nueva = []
        chl_nueva = []

        for i in range(len(fecha1) - 1):  # Añadimos el punto actual
            # print(i)
            lat_nueva.append(fecha1.iloc[i]['latitud'])
            lon_nueva.append(fecha1.iloc[i]['longitud'])
            chl_nueva.append(fecha1.iloc[i]['clorofila'])      
            if saltos.iloc[i]: # Si hay salto, añadimos un NaN
                # print('añado NAN')
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

        df_completo.loc[df_completo["clorofila"] > 1e5, "clorofila"] = np.nan

        # print(df_completo.head(5))  # vista previa
        # hacer_plot_nan(df_completo)

        # INTERPOlAR CLOROFILA
        punts_fecha1 = np.column_stack((df_completo.longitud, df_completo.latitud))
        valores_fecha1 = df_completo.clorofila.values

        # clorofila_interpolada = griddata(punts_fecha1, valores_fecha1, (lon_grid, lat_grid), method='linear')

        # print(f' hay {(~np.isnan(clorofila_interpolada)).sum()} valores no nan en la interpolacion de la fecha')
        # print(f' hay {(~np.isnan(fecha1.clorofila)).sum()} valores no nan en la  fecha real')

        # df_clorofila_mascarada = get_df_clorofila_mascarada(lon_grid, lat_grid, mascara, clorofila_interpolada, fecha)

        # if df_final.empty:
        #     print("El DataFrame está vacío.")
        #     df_final = df_clorofila_mascarada
        # else:
        #     print("El DataFrame tiene datos.")
        #     df_final = pd.concat([df_final, df_clorofila_mascarada], ignore_index=True)

        # NEAREST
        # clorofila_interpolada_nearest = griddata(punts_fecha1, valores_fecha1, (lon_grid, lat_grid), method='nearest')

        # print(f' hay {(~np.isnan(clorofila_interpolada_nearest)).sum()} valores no nan en la interpolacion de la fecha')
        # print(f' hay {(~np.isnan(fecha1.clorofila)).sum()} valores no nan en la  fecha real')

        # df_clorofila_mascarada_nearest = get_df_clorofila_mascarada(lon_grid, lat_grid, mascara, clorofila_interpolada_nearest, fecha)

        # if df_final.empty:
        #     print("El DataFrame está vacío.")
        #     df_final = df_clorofila_mascarada
        # else:
        #     print("El DataFrame tiene datos.")
        #     df_final = pd.concat([df_final, df_clorofila_mascarada], ignore_index=True)
        
        # LINEAR-NEAREST
        # clorofila_interpolada_ln = griddata(punts_fecha1, valores_fecha1, (lon_grid, lat_grid), method='linear')
        # # Dónde hay NaN
        # mask_nan = np.isnan(clorofila_interpolada_ln)
        # # Interpolar solo esos puntos usando nearest
        # clorofila_interpolada_ln[mask_nan] = griddata(punts_fecha1, valores_fecha1, (lon_grid[mask_nan], lat_grid[mask_nan]), method='nearest')

        # print(f' hay {(~np.isnan(clorofila_interpolada_ln)).sum()} valores no nan en la interpolacion de la fecha')
        # print(f' hay {(~np.isnan(fecha1.clorofila)).sum()} valores no nan en la  fecha real')

        # df_clorofila_mascarada_ln = get_df_clorofila_mascarada(lon_grid, lat_grid, mascara, clorofila_interpolada_ln, fecha)

        # if df_final.empty:
        #     print("El DataFrame está vacío.")
        #     df_final = df_clorofila_mascarada
        # else:
        #     print("El DataFrame tiene datos.")
        #     df_final = pd.concat([df_final, df_clorofila_mascarada], ignore_index=True)
        
        # MASCARA DE DISTANCIAS
        # 2. Árbol para distancias
        clorofila_interpolada_tree = griddata(punts_fecha1, valores_fecha1, (lon_grid, lat_grid), method='nearest')
        tree = cKDTree(punts_fecha1)
        distancias, _ = tree.query(np.column_stack((lon_grid.ravel(), lat_grid.ravel())))
        distancias = distancias.reshape(lon_grid.shape)

        # 3. Máscara de distancia máxima
        distancias_lat = pd.Series( [(fecha1.iloc[i+1].latitud - fecha1.iloc[i].latitud) for i in range(len(fecha1) - 1)])
        distancias_lat= distancias_lat[distancias_lat != 0]
    
        distancia_maxima = abs(distancias_lat.median()*0.75)  # Ajustar este número
        clorofila_interpolada_tree[distancias > distancia_maxima] = np.nan

        df_clorofila_mascarada_tree = get_df_clorofila_mascarada(lon_grid, lat_grid, mascara, clorofila_interpolada_tree, fecha)

        # hacer_5_plots(df1 = fecha1,
        #             df2= df_clorofila_mascarada,
        #             df3 = df_clorofila_mascarada_nearest,
        #             df4= df_clorofila_mascarada_ln,
        #             df5= df_clorofila_mascarada_tree )
        
        hacer_plot(df_clorofila_mascarada_tree)

        if df_final.empty:
            print("El DataFrame está vacío.")
            df_final = df_clorofila_mascarada_tree
        else:
            print("El DataFrame tiene datos.")
            df_final = pd.concat([df_final, df_clorofila_mascarada_tree], ignore_index=True)
    
# %%
df_final.to_csv('datos_satelite_limpios_sin10000.csv', index=False)
# %%
