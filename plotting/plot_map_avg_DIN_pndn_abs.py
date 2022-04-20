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
avg = True
mmol = True

dp_layer = 20 # choose depth layer, 25 = 50 m
#dp_avg0 = 15 # depth layers to average over 15 = 30 m
#dp_avg1 = 25 # depth layers to average over
dp_avg0 = 0 # depth layers to average over 15 = 30 m
dp_avg1 = 39 # depth layers to average over
depth = 80

# avg maps
outpath = '/data/project6/minnaho/opc_scenarios/ext_depth/'
filename = 'ext_0_80_'

savepath = './figs/maps/'

# month or season
timename = 'spring1998'

timeunit = 'days since 1997-08-01'

exp = ['PNDN_only','pndn50','pndn90']
title_exp = ['PNDN only','PNDN 50','PNDN 90']
#exp = ['l1617','FNDN_only','fndn50','fndn90']
#title_exp = ['Loads 16-17','FNDN only','FNDN 50','FNDN 90']
#exp = ['l1617']
#title_exp = ['Loads 16-17']

# roms var
var_nc = 'var'
varstr = 'DIN'
cblabel = varstr+' mmol m$^{-3}$'
if perc == True:
    cblabel = varstr+' % change'

mmolm3_to_mgl = 14./1000

# color map
#c_map = cmocean.cm.dense
#c_map = cmocean.cm.thermal
#c_map = cmocean.cm.algae
#c_map = cmocean.cm.deep
#c_map = cmocean.cm.delta
c_map = 'PRGn'
c_map1 = cmocean.cm.dense

# control scenario
cntrl_nc = outpath+filename+varstr+'_'+timename+'_PNDN_only.nc'


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
figh = 12

axis_tick_size = 16

# large pipes lat and lon
major_nc = Dataset('/data/project1/minnaho/potw_outfall_data/updated_2013_2017/major_potw_data/major_potw_1971_2017_monthly.nc','r')
lat_potw = np.array(major_nc.variables['latitude'])
lon_potw = np.array(major_nc.variables['longitude'])

coast_10m = cpf.NaturalEarthFeature('physical','coastline','10m')

# max and min of color bar
if var_nc == 'var':
    if perc == True:
        # 1998
        v_max = 15
        v_min = -15
        # 1999
        #v_max = 50
        #v_min = -30
        if avg == True:
            v_max = 50
            v_min = -50
    else:
        v_max = 1
        v_min = -1


fig,ax = plt.subplots(1,3,figsize=[figw,figh],subplot_kw=dict(projection=ccrs.PlateCarree()))

for t_i in range(len(ncfiles)):
    datanc = np.squeeze(Dataset(ncfiles[t_i],'r')['var'])
    datanc[datanc>1E10] = np.nan
    datanc[datanc<=0] = np.nan
    if avg == True:
        varrd = np.nanmean(datanc[dp_avg0:dp_avg1+1],axis=0)
    else:
        varrd = datanc[dp_layer,:,:]
    varrd[varrd>1E10] = np.nan
    varrd[varrd<=0] = np.nan
    #varplt = varrd
    
    cntrld = np.squeeze(Dataset(cntrl_nc,'r')['var'])
    cntrld[cntrld>1E10] = np.nan
    cntrld[cntrld<=0] = np.nan
    if avg == True:
        cntrlv = np.nanmean(cntrld[dp_avg0:dp_avg1+1],axis=0)
    else:
        cntrlv = cntrld[dp_layer,:,:]
    cntrlv[cntrlv>1E10] = np.nan
    cntrlv[cntrlv<=0] = np.nan

    if perc == True:    
        varplt = ((varrd - cntrlv)/cntrlv)*100
    else:    
        if 'PNDN' in ncfiles[t_i]:
            varplt = varrd
        else:
            varplt = (varrd - cntrlv)
            #varplt = (varrd - cntrlv)*mmolm3_to_mgl

    # plot maps
    #p_plot = ax.flat[t_i].pcolormesh(lon_nc,lat_nc,varplt,transform=ccrs.PlateCarree(),cmap=c_map,vmin=v_min,vmax=v_max)
    
    if 'PNDN' in ncfiles[t_i]:
        p_plot1 = ax.flat[t_i].pcolormesh(lon_nc,lat_nc,varplt,transform=ccrs.PlateCarree(),cmap=c_map1,vmin=0,vmax=5)
    else:
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
    if a_i >= 1:
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
p0 = ax.flat[0].get_position().get_points().flatten()
p1 = ax.flat[0].get_position().get_points().flatten()
cb_ax1 = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])

cb = fig.colorbar(p_plot1,cax=cb_ax1,orientation='vertical')
#cb.set_label(cblabel,fontsize=axis_tick_size)
cb.ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)

p0 = ax.flat[2].get_position().get_points().flatten()
p1 = ax.flat[2].get_position().get_points().flatten()
cb_ax = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])

cb = fig.colorbar(p_plot,cax=cb_ax,orientation='vertical')
cb.set_label(cblabel,fontsize=axis_tick_size)
cb.ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)

savename = 'map_avg_'+str(dp_avg0*2)+'_'+str(dp_avg1*2)+'_'+varstr+'_'+timename+'_'+region_name+'_'+exp[-1]+'_minus_PNDN_abs.png'

#plt.tight_layout()

fig.savefig(savepath+savename,bbox_inches='tight')
print(savename)
#plt.close()

