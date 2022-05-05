################################################
# plot maps of opc scenario
# cntrl
# compared to loads 16-17
################################################
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
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

plt.ion()

# avg maps
outpath = '/data/project3/minnaho/opc_scenarios/plotting/habitat_capacity/maps/'
filename = 'map_omega_th_1.4_'
dtstr = 'Y1999M04'
#filename = 'avg_alltime_map_omega_th_1.4'
#dtstr = ''

savepath = './figs/maps/'

#exp = ['l1617','PNDN_only','pndn50','pndn90']
#title_exp = ['Loads 16-17','PNDN only','PNDN 50','PNDN 90']

exp = ['l1617']
title_exp = ['Loads 16-17'] 

# roms var
cblabel = '% change in habitat capacity'

# color map
#c_map = cmocean.cm.delta
#c_map = 'PRGn'
c_map = cmocean.cm.balance_r

# outputs
ncfiles = []
for e_i in range(len(exp)):
    ncfiles.append(outpath+filename+dtstr+'_'+exp[e_i]+'_cntrl.nc')

# grid path
grid_nc = l2grid.grid_nc
lat_nc = l2grid.lat_nc
lon_nc = l2grid.lon_nc
h_nc = l2grid.h_nc

# LA/OC region
#region_name = 'laoc'
#lat_min = 33.5
#lat_max = 34.1
#lon_min = -118.9
#lon_max = -117.82

# full grid
region_name = 'grid'
lat_min = np.nanmin(lat_nc)
lat_max = np.nanmax(lat_nc)
lon_min = np.nanmin(lon_nc)
lon_max = np.nanmax(lon_nc)

extent = [lon_min,lon_max,lat_min,lat_max]

figw = 12
figh = 10

axis_tick_size = 16

# large pipes lat and lon
major_nc = Dataset('/data/project1/minnaho/potw_outfall_data/updated_2013_2017/major_potw_data/major_potw_1971_2017_monthly.nc','r')
lat_potw = np.array(major_nc.variables['latitude'])
lon_potw = np.array(major_nc.variables['longitude'])

coast_10m = cpf.NaturalEarthFeature('physical','coastline','10m')

# max and min of color bar
v_max = 40
v_min = -40
#v_max = 20
#v_min = -20
#v_max = 5
#v_min = -5

t_i = 0

fig,ax = plt.subplots(1,1,figsize=[figw,figh],subplot_kw=dict(projection=ccrs.PlateCarree()))

datanc = Dataset(ncfiles[t_i],'r')
varplt = np.squeeze(np.array(datanc['habitat_cap']))

# plot maps
p_plot = ax.pcolormesh(lon_nc,lat_nc,varplt,transform=ccrs.PlateCarree(),cmap=c_map,norm=mcolors.DivergingNorm(0),vmin=v_min,vmax=v_max)
#p_plot = ax.pcolormesh(lon_nc,lat_nc,varplt,transform=ccrs.PlateCarree(),cmap=c_map,norm=mcolors.DivergingNorm(0))

# contour
c_plotneg = ax.contour(lon_nc,lat_nc,varplt,[-10],colors='red')
c_plotpos = ax.contour(lon_nc,lat_nc,varplt,[10],colors='blue')

ax.set_title(title_exp[t_i],fontsize=axis_tick_size)

# mark pipe location
for l_i in range(len(lon_potw)):
    ax.scatter(lon_potw[l_i],lat_potw[l_i],marker='o',facecolors='none',edgecolors='gray',s=100)
# other grid stuff
ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.add_feature(coast_10m,facecolor='None',edgecolor='k')
ax.add_feature(cpf.BORDERS,facecolor='None',edgecolor='k')
ax.set_extent(extent)
# lat/lon axes
gl = ax.gridlines(draw_labels=True,linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlabel_style = {'size':axis_tick_size}
gl.ylabel_style = {'size':axis_tick_size}
#step_lon = .4
#step_lat = .2
step_lon = 2
step_lat = 1
gl.xlocator = mticker.FixedLocator(np.arange(lon_min-step_lon,lon_max+step_lon,step_lon).astype(int))
gl.ylocator = mticker.FixedLocator(np.arange(lat_min-step_lat,lat_max+step_lat,step_lat).astype(int))
#gl.ylocator = mticker.FixedLocator([33.4, 33.5, 33.6, 33.7, 33.8, 33.9, 34. , 34.1, 34.2])
#gl.ylocator = mticker.FixedLocator([33.5, 33.7, 33.9, 34.1])
gl.yformatter = LATITUDE_FORMATTER
gl.xformatter = LONGITUDE_FORMATTER

# colorbar
p0 = ax.get_position().get_points().flatten()
p1 = ax.get_position().get_points().flatten()
cb_ax = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])

cb = fig.colorbar(p_plot,cax=cb_ax,orientation='vertical',format='%.1i')
cb.set_label(cblabel,fontsize=axis_tick_size)
cb.ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)

savename = filename+dtstr+'_'+exp[-1]+'_cntrl.png'

plt.tight_layout()

fig.savefig(savepath+savename,bbox_inches='tight')
print(savename)
#plt.close()

