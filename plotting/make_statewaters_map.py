import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import cmocean as cmocean
import cartopy.crs as ccrs
import cartopy.feature as cpf
from cartopy.io import shapereader as shpreader
import geopandas
from osgeo import ogr

fig_w = 15
fig_h = 12

plt.ion()

lat_min = 32.4
lat_max = 34.6
lon_min = -120.5
lon_max = -117

extent = [lon_min,lon_max,lat_min,lat_max]

coast_10m = cpf.NaturalEarthFeature('physical','coastline','10m')

fig,ax = plt.subplots(1,1,figsize=[fig_w,fig_h],subplot_kw=dict(projection=ccrs.PlateCarree()))
ax.add_feature(coast_10m,facecolor='None',edgecolor='k')
ax.set_extent(extent)



shpfile = '../ca_statewaters_shapefile/'
shp = geopandas.read_file(shpfile+'CA_cst3nm.shp')
border = shp.geometry[0]

cashp = cpf.ShapelyFeature(shpreader.Reader(shpfile).geometries(),ccrs.PlateCarree(),edgecolor='dodgerblue',facecolor='None',alpha=0.5)

ax.add_feature(border)

ds = ogr.Open(shpfile+'CA_cst3nm.shp')
lyr = ds.GetLayerByName('CA_cst3nm')
env = lyr[0].GetGeometryRef().GetEnvelope()

#https://gis.stackexchange.com/questions/3821/how-can-i-convert-a-shapefile-to-lat-and-lon-boundaries


