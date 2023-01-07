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

# avg maps
outpath = '/data/project6/minnaho/opc_scenarios/bgc_flux/'
filename = 'concat_fullts_int_avg_100m_50m_'

savepath = './figs/maps/'

#exp = ['PNDN_only','pndn50','pndn90']
#title_exp = ['PNDN only','PNDN 50','PNDN 90']
#exp = ['l1617','FNDN_only','fndn50','fndn90']
#title_exp = ['Loads 16-17','FNDN only','FNDN 50','FNDN 90']
#exp = ['l1617']
#title_exp = ['Loads 16-17']
#exp = ['cntrl_initap_realistic','fulll_2012_2017','PNDN_only_realistic','pndn50_fixriver','pndn90_fixriver','FNDN_only_realistic','fndn50_fixriver','fndn90_fixriver']
exp = ['cntrl_initap_realistic','fulll_2012_2017']
title_exp = ['CTRL','ANTH']

# roms var
var_nc = 'var_int'
varstr = 'npp'
cblabel = 'NPP mmol m$^{-2}$ d$^{-1}$'
if perc == True:
    cblabel = varstr+' % change'

s2d = 86400

# color map
#c_map = cmocean.cm.dense
#c_map = cmocean.cm.thermal
#c_map = cmocean.cm.algae
#c_map = cmocean.cm.deep
c_map = cmocean.cm.delta
#c_map = 'PRGn'
#c_map = cmocean.cm.balance
#c_map1 = cmocean.cm.algae

# control scenario
cntrl_nc = outpath+filename+'fulll_2012_2017'+'_'+varstr+'.nc'
# subtract from cntrl
cntrld = np.squeeze(Dataset(cntrl_nc,'r')[var_nc])
cntrld[cntrld>1E10] = np.nan
cntrld[cntrld<=0] = np.nan
# time average
cntrlavg = np.nanmean(cntrld,axis=0)

# outputs
ncfiles = []
for e_i in range(len(exp)):
    ncfiles.append(outpath+filename+exp[e_i]+'_'+varstr+'.nc')

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

region_name = 'wider'
lat_min = 31.9
lat_max = 34.6
lon_min = -120.7
lon_max = -117

# full grid
#region_name = 'grid'
#lat_min = np.nanmin(lat_nc)
#lat_max = np.nanmax(lat_nc)
#lon_min = np.nanmin(lon_nc)
#lon_max = np.nanmax(lon_nc)

extent = [lon_min,lon_max,lat_min,lat_max]

figw = 14
figh = 7

axfont = 16

# large pipes lat and lon
major_nc = Dataset('/data/project1/minnaho/potw_outfall_data/updated_2013_2017/major_potw_data/major_potw_1971_2017_monthly.nc','r')
lat_potw = np.array(major_nc.variables['latitude'])
lon_potw = np.array(major_nc.variables['longitude'])

coast_10m = cpf.NaturalEarthFeature('physical','coastline','10m')

# max and min of color bar
v_max = 150
v_min = 50


fig,ax = plt.subplots(1,2,figsize=[figw,figh],subplot_kw=dict(projection=ccrs.PlateCarree()))

for e_i in range(len(exp)):
    print(title_exp[e_i])
    # convert s^-1 to d^-1
    datanc = np.squeeze(Dataset(ncfiles[e_i],'r')[var_nc])*86400
    datanc[datanc>1E10] = np.nan
    datanc[datanc<=0] = np.nan
    # time average
    # full 2 years
    dataavg = np.nanmean(datanc,axis=0)
    # spring
    #dataavg = np.nanmean((np.concatenate((datanc[5:8],datanc[17:20]))),axis=0)
    # summer
    #dataavg = np.nanmean((np.concatenate((datanc[8:11],datanc[20:23]))),axis=0)
    # summer
    #dataavg = np.nanmean((np.concatenate((datanc[2:5],datanc[14:17]))),axis=0)

    #datad = dataavg - cntrlavg
    #varplt = datad
    varplt = dataavg
    #p_plot = ax.flat[e_i].pcolormesh(lon_nc,lat_nc,varplt,transform=ccrs.PlateCarree(),cmap=c_map,norm=mcolors.DivergingNorm(0))
    #p_plot = ax.flat[e_i].pcolormesh(lon_nc,lat_nc,varplt,transform=ccrs.PlateCarree(),vmin=v_min,vmax=v_max,cmap=c_map,norm=mcolors.DivergingNorm(0))
    p_plot = ax.flat[e_i].pcolormesh(lon_nc,lat_nc,varplt,transform=ccrs.PlateCarree(),cmap=c_map,norm=mcolors.DivergingNorm(vmin=0,vcenter=80,vmax=v_max))
    
    ax.flat[e_i].set_title(title_exp[e_i],fontsize=axfont)

for a_i in range(len(ax.flat)):
    # mark pipe location
    for l_i in range(len(lon_potw)):
        ax.flat[a_i].scatter(lon_potw[l_i],lat_potw[l_i],marker='o',facecolors='none',edgecolors='blue',s=100)
    # other grid stuff
    ax.flat[a_i].add_feature(coast_10m,facecolor='None',edgecolor='k')
    ax.flat[a_i].add_feature(cpf.BORDERS,facecolor='None',edgecolor='k')
    ax.flat[a_i].set_extent(extent)
    gl = ax.flat[a_i].gridlines(draw_labels={'bottom':'x','left':'y'},xlabel_style=dict(size=axfont),ylabel_style=dict(size=axfont))
    if a_i >= 1:
        gl.left_labels=False

    #ax.flat[a_i].tick_params(axis='both',which='major',labelsize=axfont)

# colorbar
p0 = ax.flat[1].get_position().get_points().flatten()
p1 = ax.flat[1].get_position().get_points().flatten()
cb_ax = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])

cb = fig.colorbar(p_plot,cax=cb_ax,orientation='vertical')
cb.set_label(cblabel,fontsize=axfont)
cb.ax.tick_params(axis='both',which='major',labelsize=axfont)


savename = 'map_2years_cntrl_anth_100m_int_'+varstr+'_'+region_name+'_abs_fullts.png'

#plt.tight_layout()

fig.savefig(savepath+savename,bbox_inches='tight')
print(savename)

