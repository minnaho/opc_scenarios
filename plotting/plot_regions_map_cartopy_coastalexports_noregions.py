import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import cmocean as cmocean
import h5py
import cartopy.crs as ccrs
import cartopy.feature as cpf
from cartopy.io import shapereader as shpreader
import cartopy.io.img_tiles as cimgt

plt.ion()

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


# grid
grid_path = '/data/project5/kesf/ROMS/L2_SCB/roms_grd.nc'
grid_nc = Dataset(grid_path,'r')

lat_nc = np.array(grid_nc.variables['lat_rho'])
lon_nc = np.array(grid_nc.variables['lon_rho'])

lat_min = 32.4
lat_max = 34.6
lon_min = -120.5
lon_max = -117

# plot
axis_tick_size = 16
# latitudes to draw
parallels = np.arange(0,90,1)
# longitudes to draw
meridians = np.arange(180,360,1)

extent = [lon_min,lon_max,lat_min,lat_max]
rivers_10m = cpf.NaturalEarthFeature('physical','rivers_lake_centerlines','10m')
coast_10m = cpf.NaturalEarthFeature('physical','coastline','10m')
wsheds = shpreader.Reader(shp_path+'basin_arcgis/wribasin.shp')

fig_w = 15
fig_h = 12

fig,ax = plt.subplots(1,1,figsize=[fig_w,fig_h],subplot_kw=dict(projection=ccrs.PlateCarree()))
ax.add_feature(coast_10m,facecolor='None',edgecolor='k')
ax.add_feature(cpf.BORDERS,facecolor='None',edgecolor='k')
ax.set_extent(extent)
gl = ax.gridlines(draw_labels={'bottom':'x','left':'y'},xlabel_style=dict(size=axfont),ylabel_style=dict(size=axfont))

#ax.set_extent(extent)
#gl = ax.gridlines(draw_labels=True,linestyle='--')
#gl.xlabels_top = False
#gl.ylabels_right = False
#gl.xlabel_style = {'size':20}
#gl.ylabel_style = {'size':20}
#ax.add_feature(coast_10m,facecolor='None',edgecolor='k')

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

fig.savefig('./figs/2years/inputs_map.png',bbox_inches='tight')
