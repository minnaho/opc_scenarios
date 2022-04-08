################################################
# plot maps of opc scenarios - 8 maps total
# 1 month or season
# 7 OPC scenarios
# 1 CTRL scenario
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

perc = False
avg = False

dp_layer = 40 # choose depth layer, 25 = 50 m
depth = 80

# avg maps
outpath = '/data/project6/minnaho/opc_scenarios/ext_depth_200/'
filename = 'ext_0_200_'

savepath = './figs/maps/'

# month or season
timename = 'spring1998'

timeunit = 'days since 1997-08-01'

#exp = ['l1617','PNDN_only','pndn50','pndn90']
#title_exp = ['Loads 16-17','PNDN only','PNDN 50','PNDN 90']
exp = ['PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
title_exp = ['PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']
#exp = ['l1617']
#title_exp = ['Loads 16-17']

# roms var
var_nc = 'omega'
varstr = 'omega_pH'
cblabel = var_nc+' change'
if perc == True:
    cblabel = var_nc+' % change'

# color map
#c_map = cmocean.cm.dense
#c_map = cmocean.cm.thermal
#c_map = cmocean.cm.algae
#c_map = cmocean.cm.deep
#c_map = cmocean.cm.delta
c_map = 'PRGn'

# control scenario
cntrl_nc = outpath+filename+varstr+'_'+timename+'_l1617.nc'


# outputs
ncfiles = []
for e_i in range(len(exp)):
    ncfiles.append(outpath+filename+varstr+'_'+timename+'_'+exp[e_i]+'.nc')

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

figw = 18
figh = 10

axis_tick_size = 16

# large pipes lat and lon
major_nc = Dataset('/data/project1/minnaho/potw_outfall_data/updated_2013_2017/major_potw_data/major_potw_1971_2017_monthly.nc','r')
lat_potw = np.array(major_nc.variables['latitude'])
lon_potw = np.array(major_nc.variables['longitude'])

coast_10m = cpf.NaturalEarthFeature('physical','coastline','10m')

# max and min of color bar
if perc == True:
    # 1998
    v_max = 10
    v_min = -10
    # 1999
    #v_max = 50
    #v_min = -30
    if avg == True:
        v_max = 15
        v_min = -15
else:
    if var_nc == 'omega':
        v_max = 0.1
        v_min = -0.1
    if var_nc == 'pH':
        v_max = 0.02
        v_min = -0.02


fig,ax = plt.subplots(2,3,figsize=[figw,figh],subplot_kw=dict(projection=ccrs.PlateCarree()))

for t_i in range(len(ncfiles)):
    datanc = Dataset(ncfiles[t_i],'r')
    if avg == True:
        varrd = np.nanmean(np.squeeze(np.array(datanc[var_nc])),axis=0)
    else:
        varrd = np.squeeze(np.array(datanc[var_nc]))[dp_layer,:,:]
    varrd[varrd>1E10] = np.nan
    varrd[varrd<=0] = np.nan
    
    cntrld = Dataset(cntrl_nc,'r')
    if avg == True:
        cntrlv = np.nanmean(np.squeeze(np.array(cntrld[var_nc])),axis=0)
    else:
        cntrlv = np.squeeze(np.array(cntrld[var_nc]))[dp_layer,:,:]
    cntrlv[cntrlv>1E10] = np.nan
    cntrlv[cntrlv<=0] = np.nan
    
    if perc == True:    
        varplt = ((varrd - cntrlv)/cntrlv)*100
    else:    
        varplt = varrd - cntrlv

    # plot maps
    p_plot = ax.flat[t_i].pcolormesh(lon_nc,lat_nc,varplt,transform=ccrs.PlateCarree(),cmap=c_map,norm=mcolors.DivergingNorm(0),vmin=v_min,vmax=v_max)
    #p_plot = ax.flat[t_i].pcolormesh(lon_nc,lat_nc,varplt,transform=ccrs.PlateCarree(),cmap=c_map,norm=mcolors.DivergingNorm(0))
    
    #fig.suptitle(timename,fontsize=axis_tick_size)
    
    ax.flat[t_i].set_title(title_exp[t_i],fontsize=axis_tick_size)
    
    
for a_i in range(len(ax.flat)):
    # mark pipe location
    for l_i in range(len(lon_potw)):
        ax.flat[a_i].scatter(lon_potw[l_i],lat_potw[l_i],marker='o',facecolors='none',edgecolors='blue',s=100)
    # other grid stuff
    ax.flat[a_i].tick_params(axis='both',which='major',labelsize=axis_tick_size)
    ax.flat[a_i].yaxis.set_ticks_position('both')
    ax.flat[a_i].xaxis.set_ticks_position('both')
    ax.flat[a_i].add_feature(coast_10m,facecolor='None',edgecolor='k')
    ax.flat[a_i].add_feature(cpf.BORDERS,facecolor='None',edgecolor='k')
    ax.flat[a_i].set_extent(extent)
    # lat/lon axes
    gl = ax.flat[a_i].gridlines(draw_labels=True,linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlabel_style = {'size':axis_tick_size}
    gl.ylabel_style = {'size':axis_tick_size}
    if a_i >= 1 and a_i != 3:
        gl.ylabels_left = False
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
p0 = ax.flat[2].get_position().get_points().flatten()
p1 = ax.flat[5].get_position().get_points().flatten()
cb_ax = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])

#cb = fig.colorbar(p_plot,cax=cb_ax,orientation='vertical',format='%.1f')
cb = fig.colorbar(p_plot,cax=cb_ax,orientation='vertical')
cb.set_label(cblabel,fontsize=axis_tick_size)
cb.ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)

if perc == True:
    savename = 'map_'+str(depth)+'m_'+var_nc+'_'+timename+'_'+region_name+'_'+exp[-1]+'_l1617_perc.png'
    if avg == True:
        savename = 'map_avg_0_200_'+var_nc+'_'+timename+'_'+region_name+'_'+exp[-1]+'_l1617_perc.png'
else:
    savename = 'map_'+str(depth)+'m_'+var_nc+'_'+timename+'_'+region_name+'_'+exp[-1]+'_l1617_abs.png'

plt.tight_layout()

fig.savefig(savepath+savename,bbox_inches='tight')
print(savename)
#plt.close()

