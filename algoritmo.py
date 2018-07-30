#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#######################################
# Script que permite la generación
# automática de mapas de particulas en el
# aire
# Author: Jorge Mauricio
# Email: jorge.ernesto.mauricio@gmail.com
# Date: Created on Thu Sep 28 08:38:15 2017
# Version: 1.0
#######################################
"""

# librerías
import os
import pandas as pd
import numpy as np
import h5py
import requests
from bs4 import BeautifulSoup
import urllib.request
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import time
from time import gmtime, strftime

# límites Lat y Long
LONG_MIN = -115
LONG_MAX = -111
LAT_MIN = 29
LAT_MAX = 32

PATH = "/home/jorge/Documents/Research/imagenesCaborca"

array_URLs = ["https://acdisc.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level3/OMNO2d.003/2018/",
              "https://acdisc.gsfc.nasa.gov/data/Aura_OMI_Level3/OMDOAO3e.003/2018/",
              "https://acdisc.gsfc.nasa.gov/data/Aura_OMI_Level3/OMSO2e.003/2018/",
              "https://acdisc.gsfc.nasa.gov/data/Aura_OMI_Level3/OMTO3e.003/2018/",
              "https://acdisc.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level3/OMAEROe.003/2018/"]

#array_URLs = ["https://acdisc.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level3/OMNO2d.003/2018/"]

array_Archivo = []



# función main
def main():
    # fecha de la descarga
    # descarga de información
    descarga_de_archivos()
    procesamientoNO2()
    procesamientoO3()
    procesamientoSO2()
    procesamientoTO3()
    procesamientoAERO()


# función para procesar NO2
def procesamientoNO2():
    # clear plt
    plt.clf()
    # Open file.
    FILE_NAME = array_Archivo[0]
    DATAFIELD_NAME = 'HDFEOS/GRIDS/ColumnAmountNO2/Data Fields/ColumnAmountNO2'
    with h5py.File(FILE_NAME, mode='r') as f:
        # Read dataset.
        dset = f[DATAFIELD_NAME]
        data = dset[:]

        # Handle fill value.
        data[data == dset.fillvalue] = np.nan
        data = np.ma.masked_where(np.isnan(data), data)

        # Get attributes needed for the plot.
        # String attributes actually come in as the bytes type and should
        # be decoded to UTF-8 (python3).
        title = dset.attrs['Title'].decode()
        units = dset.attrs['Units'].decode()

        # There is no geolocation data, so construct it ourselves.
    longitude = np.arange(0., 1440.0) * 0.25 - 180 + 0.125
    latitude = np.arange(0., 720.0) * 0.25 - 90 + 0.125


    # leer coordenadas
    dataEstaciones = pd.read_csv("{}/data/coordenadas_estaciones.csv".format(PATH))
    xC = np.array(dataEstaciones['Long'])
    yC = np.array(dataEstaciones['Lat'])

    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=LAT_MIN, urcrnrlat = LAT_MAX,
                llcrnrlon=LONG_MIN, urcrnrlon = LONG_MAX)

    m.scatter(xC, yC, latlon=True, s=1, marker='o', color='r', zorder=25)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180., 45.), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, data, latlon=True, cmap='jet')

    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, title))

    fig = plt.gcf()
    # plt.show()

    pngfile = "{}.png".format(basename)
    fig.savefig(pngfile, dpi=600)

# función para procesar 03
def procesamientoO3():
    # clear plt
    plt.clf()
    # Open file.
    FILE_NAME = array_Archivo[1]
    DATAFIELD_NAME = 'HDFEOS/GRIDS/ColumnAmountO3/Data Fields/ColumnAmountO3'
    with h5py.File(FILE_NAME, mode='r') as f:
        # Read dataset.
        dset = f[DATAFIELD_NAME]
        data = dset[:]

        # Handle fill value.
        data[data == dset.fillvalue] = np.nan
        data = np.ma.masked_where(np.isnan(data), data)

        # Get attributes needed for the plot.
        # String attributes actually come in as the bytes type and should
        # be decoded to UTF-8 (python3).
        title = dset.attrs['Title'].decode()
        units = dset.attrs['Units'].decode()

        # There is no geolocation data, so construct it ourselves.
    longitude = np.arange(0., 1440.0) * 0.25 - 180 + 0.125
    latitude = np.arange(0., 720.0) * 0.25 - 90 + 0.125


    # leer coordenadas
    dataEstaciones = pd.read_csv("{}/data/coordenadas_estaciones.csv".format(PATH))
    xC = np.array(dataEstaciones['Long'])
    yC = np.array(dataEstaciones['Lat'])

    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=LAT_MIN, urcrnrlat = LAT_MAX,
                llcrnrlon=LONG_MIN, urcrnrlon = LONG_MAX)

    m.scatter(xC, yC, latlon=True, s=1, marker='o', color='r', zorder=25)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180., 45.), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, data, latlon=True, cmap='jet')

    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, title))

    fig = plt.gcf()
    # plt.show()

    pngfile = "{}.png".format(basename)
    fig.savefig(pngfile, dpi=600)

# función para procesar SO2
def procesamientoSO2():
    # clear plt
    plt.clf()
    # Open file.
    FILE_NAME = array_Archivo[2]
    DATAFIELD_NAME = 'HDFEOS/GRIDS/OMI Total Column Amount SO2/Data Fields/ColumnAmountSO2_PBL'
    with h5py.File(FILE_NAME, mode='r') as f:
        # Read dataset.
        dset = f[DATAFIELD_NAME]
        data = dset[:]

        # Handle fill value.
        data[data == dset.fillvalue] = np.nan
        data = np.ma.masked_where(np.isnan(data), data)

        # Get attributes needed for the plot.
        # String attributes actually come in as the bytes type and should
        # be decoded to UTF-8 (python3).
        title = dset.attrs['Title'].decode()
        units = dset.attrs['Units'].decode()

        # There is no geolocation data, so construct it ourselves.
    longitude = np.arange(0., 1440.0) * 0.25 - 180 + 0.125
    latitude = np.arange(0., 720.0) * 0.25 - 90 + 0.125


    # leer coordenadas
    dataEstaciones = pd.read_csv("{}/data/coordenadas_estaciones.csv".format(PATH))
    xC = np.array(dataEstaciones['Long'])
    yC = np.array(dataEstaciones['Lat'])

    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=LAT_MIN, urcrnrlat = LAT_MAX,
                llcrnrlon=LONG_MIN, urcrnrlon = LONG_MAX)

    m.scatter(xC, yC, latlon=True, s=1, marker='o', color='r', zorder=25)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180., 45.), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, data, latlon=True, cmap='jet')

    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, title))

    fig = plt.gcf()
    # plt.show()

    pngfile = "{}.png".format(basename)
    fig.savefig(pngfile, dpi=600)

# función para procesar TO3
def procesamientoTO3():
    # clear plt
    plt.clf()
    # Open file.
    FILE_NAME = array_Archivo[3]
    DATAFIELD_NAME = 'HDFEOS/GRIDS/OMI Column Amount O3/Data Fields/ColumnAmountO3'
    with h5py.File(FILE_NAME, mode='r') as f:
        # Read dataset.
        dset = f[DATAFIELD_NAME]
        data = dset[:]

        # Handle fill value.
        data[data == dset.fillvalue] = np.nan
        data = np.ma.masked_where(np.isnan(data), data)

        # Get attributes needed for the plot.
        # String attributes actually come in as the bytes type and should
        # be decoded to UTF-8 (python3).
        title = dset.attrs['Title'].decode()
        units = dset.attrs['Units'].decode()

        # There is no geolocation data, so construct it ourselves.
    longitude = np.arange(0., 1440.0) * 0.25 - 180 + 0.125
    latitude = np.arange(0., 720.0) * 0.25 - 90 + 0.125


    # leer coordenadas
    dataEstaciones = pd.read_csv("{}/data/coordenadas_estaciones.csv".format(PATH))
    xC = np.array(dataEstaciones['Long'])
    yC = np.array(dataEstaciones['Lat'])

    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=LAT_MIN, urcrnrlat = LAT_MAX,
                llcrnrlon=LONG_MIN, urcrnrlon = LONG_MAX)

    m.scatter(xC, yC, latlon=True, s=1, marker='o', color='r', zorder=25)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180., 45.), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, data, latlon=True, cmap='jet')

    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, title))

    fig = plt.gcf()
    # plt.show()

    pngfile = "{}.png".format(basename)
    fig.savefig(pngfile, dpi=600)

# función para procesar TO3
def procesamientoAERO():
    # clear plt
    plt.clf()
    # Open file.
    FILE_NAME = array_Archivo[4]
    DATAFIELD_NAME = 'HDFEOS/GRIDS/ColumnAmountAerosol/Data Fields/UVAerosolIndex'
    with h5py.File(FILE_NAME, mode='r') as f:
        # Read dataset.
        dset = f[DATAFIELD_NAME]
        data = dset[:]

        # Handle fill value.
        data[data == dset.fillvalue] = np.nan
        data = np.ma.masked_where(np.isnan(data), data)

        # Get attributes needed for the plot.
        # String attributes actually come in as the bytes type and should
        # be decoded to UTF-8 (python3).
        title = dset.attrs['Title'].decode()
        units = dset.attrs['Units'].decode()

        # There is no geolocation data, so construct it ourselves.
    longitude = np.arange(0., 1440.0) * 0.25 - 180 + 0.125
    latitude = np.arange(0., 720.0) * 0.25 - 90 + 0.125


    # leer coordenadas
    dataEstaciones = pd.read_csv("{}/data/coordenadas_estaciones.csv".format(PATH))
    xC = np.array(dataEstaciones['Long'])
    yC = np.array(dataEstaciones['Lat'])

    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat=LAT_MIN, urcrnrlat = LAT_MAX,
                llcrnrlon=LONG_MIN, urcrnrlon = LONG_MAX)

    m.scatter(xC, yC, latlon=True, s=1, marker='o', color='r', zorder=25)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 180., 45.), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, data, latlon=True, cmap='jet')

    cb = m.colorbar()
    cb.set_label(units)

    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n{1}'.format(basename, title))

    fig = plt.gcf()
    # plt.show()

    pngfile = "{}.png".format(basename)
    fig.savefig(pngfile, dpi=600)

# función descarga de archivos
def descarga_de_archivos():
    # fecha de la descarga
    fechaPronostico = strftime("%Y-%m-%d")

    # cambiar a carpeta data
    os.chdir("data")

    # crear directorio de fecha de descarga
    os.mkdir(fechaPronostico)

    # cambiar de directorio a data
    os.chdir(fechaPronostico)

    # ciclo para la descarga de información
    for URL in array_URLs:
        # generar la consulta de información
        r = requests.get(URL)
        # parsear el html para determinar los links a descargar
        soup = BeautifulSoup(r.text, "html.parser")
        # crear un array para guardar los links
        array_links = []

        # ciclo para filtrar los links con información
        for link in soup.find_all("a"):
            array_links.append(link.get("href"))

        # nombre del archivo a descargar
        nombre_archivo = array_links[-5]

        # imprimir el nombre del archivo
        print(nombre_archivo)

        # guardar el nombre del archivo para el post procesamiento
        array_Archivo.append(nombre_archivo)

        # generar el url para la descarga de información
        URL_DESCARGA = "{}{}".format(URL, nombre_archivo)

        # print url de descarga
        print(URL_DESCARGA)

        os.system("wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies {}".format(URL_DESCARGA))

if __name__ == '__main__':
    main()
