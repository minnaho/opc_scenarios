# plot domain 300 m model
import sys
import os
sys.path.append('/data/project3/minnaho/global/')
import l2grid
import numpy as np
from netCDF4 import Dataset,num2date
import glob as glob
import ROMS_depths as depths
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import cmocean
import datetime as datetime
import calendar
import cartopy.crs as ccrs
import cartopy.feature as cpf
import cartopy.io.img_tiles as cimgt
from cartopy.io import shapereader as shpreader


c_map = cmocean.cm.deep

stamen_terrain = cimgt.Stamen('terrain-background')

shp_path = '/data/project1/minnaho/inputs_update_plotting/paper_data_inputs/'

# inputs
major_nc = Dataset('/data/project1/minnaho/potw_outfall_data/updated_2013_2017/major_potw_data/major_potw_1971_2017_monthly.nc','r')
lats_major_potw = np.array(major_nc.variables['latitude'])
lons_major_potw = np.array(major_nc.variables['longitude'])

minor_nc = Dataset('/data/project1/minnaho/potw_outfall_data/updated_2013_2017/minor_potw_data/minor_potw_1997_2017_monthly.nc','r')
lats_minor_potw = np.array(minor_nc.variables['latitude'])
lons_minor_potw = np.array(minor_nc.variables['longitude'])

river_nc = Dataset('/data/project1/minnaho/river_data/updated_2013_2017/rivers_1997_2017_monthly.nc','r')

lats_river = np.array(river_nc.variables['latitude'])
lons_river = np.array(river_nc.variables['longitude'])


# grid path
grid_nc = l2grid.grid_nc
lat_nc = l2grid.lat_nc
lon_nc = l2grid.lon_nc
h_nc = l2grid.h_nc
pm_nc = l2grid.pm_nc
pn_nc = l2grid.pn_nc
mask_nc = l2grid.mask_nc

coast_10m = cpf.NaturalEarthFeature('physical','coastline','10m')

# full grid
lat_min = np.nanmin(lat_nc)
lat_max = np.nanmax(lat_nc)
lon_min = np.nanmin(lon_nc)
lon_max = np.nanmax(lon_nc)

# mask 15 km to contour
region_mask = Dataset('/data/project1/minnaho/make_masks/mask_scb.nc','r')
mask_cst = np.array(region_mask.variables['mask_coast'])
#mask_cst[mask_cst==0] = np.nan
mask_cst[np.isnan(mask_cst)] = 0
mask_onshore = mask_cst


extent = [lon_min,lon_max,lat_min,lat_max]

axfont = 16

figw = 10
figh = 8

h_nc[h_nc==3] = np.nan

plt.ion()

fig,ax = plt.subplots(1,1,figsize=[figw,figh],subplot_kw=dict(projection=ccrs.PlateCarree()))
#p_plot = ax.pcolormesh(lon_nc,lat_nc,h_nc,transform=ccrs.PlateCarree(),cmap=c_map)
# box around domain
ax.plot(lon_nc[0],lat_nc[0],color='k')
ax.plot(lon_nc[-1],lat_nc[-1],color='k')
ax.plot(lon_nc[:,-1],lat_nc[:,-1],color='k')
ax.plot(lon_nc[:,0],lat_nc[:,0],color='k')
ax.contour(lon_nc,lat_nc,mask_onshore,colors='lightgray')

ax.add_feature(coast_10m,facecolor='None',edgecolor='k')
ax.add_feature(cpf.BORDERS,facecolor='None',edgecolor='k')
ax.set_extent(extent)
gl = ax.gridlines(draw_labels={'bottom':'x','left':'y'},xlabel_style=dict(size=axfont),ylabel_style=dict(size=axfont))

# major rivers
shpfile = shp_path+'MajorRiversAndCreeks.shp'
rivershp = cpf.ShapelyFeature(shpreader.Reader(shpfile).geometries(),ccrs.PlateCarree(),edgecolor='dodgerblue',facecolor='None',alpha=0.5)
ax.add_feature(rivershp)

# streams
shpfile = shp_path+'fromabel_stream/Streams.shp'
rivershp = cpf.ShapelyFeature(shpreader.Reader(shpfile).geometries(),ccrs.PlateCarree(),edgecolor='dodgerblue',facecolor='None',alpha=0.5)
ax.add_feature(rivershp)

# river dots
m_size = 40
ax.scatter(lons_river,lats_river,s=m_size,marker='^',facecolors='lightgreen',edgecolor='green',lw=1,label='River',zorder=10)

# POTWs
m_size = 150
maj_potw_plt = ax.scatter(lons_major_potw,lats_major_potw,s=m_size,marker='o',facecolors='none',edgecolor='gold',lw=3,label='Large POTW')
#maj_potw_plt = ax.scatter(lons_major_potw,lats_major_potw,s=m_size,marker='o',facecolors='none',edgecolor='blue',lw=3)
min_potw_plt = ax.scatter(lons_minor_potw,lats_minor_potw,s=m_size,marker='s',facecolors='none',edgecolor='k',lw=2,label='Small POTW')

leg_size = 16
ax.legend(loc='lower left',fontsize=leg_size,labelspacing=1)

#p0 = ax.get_position().get_points().flatten()
#p1 = ax.get_position().get_points().flatten()
#cb_ax = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])
#cb = fig.colorbar(p_plot,cax=cb_ax,orientation='vertical')
#cb.set_label('Bathymetry (m)',fontsize=axfont)
#cb.ax.tick_params(axis='both',which='major',labelsize=axfont)
fig.savefig('./figs/2years/domain.png',bbox_inches='tight')
