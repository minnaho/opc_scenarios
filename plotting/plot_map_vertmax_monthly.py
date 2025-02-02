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

#plt.ion()

perc = False

outpath = '/data/project6/ROMS/L2SCB_OPC/'
filename = 'l2_scb_avg.'

savepath = './figs/maps/'

timeunit = 'days since 1997-08-01'

#exp = ['PNDN_only','pndn50','pndn90']
#title_exp = ['PNDN only','PNDN 50','PNDN 90']
#exp = ['l1617','FNDN_only','fndn50','fndn90']
#title_exp = ['Loads 16-17','FNDN only','FNDN 50','FNDN 90']
#exp = ['l1617']
#title_exp = ['Loads 16-17']
exp = ['loads1617','PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
title_exp = ['Loads 16-17','PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']

# choose years
start_year = 1997
end_year = 1999

# choose months between 1 and 12
start_month = 11
end_month = 11


# roms var
var_nc = 'NH4'
varstr = 'NH4'
cblabel = varstr+' mmol m$^{-3}$'
if perc == True:
    cblabel = varstr+' % change'

s2d = 86400

# color map
c_map = cmocean.cm.dense
#c_map = cmocean.cm.thermal
#c_map = cmocean.cm.algae
#c_map = cmocean.cm.deep
#c_map = cmocean.cm.delta
#c_map = 'PRGn'
c_map1 = cmocean.cm.ice_r

# control scenario
cntrlpath = '/data/project3/minnaho/postprocessing/cntrl_monthly/'

# grid path
grid_nc = l2grid.grid_nc
lat_nc = l2grid.lat_nc
lon_nc = l2grid.lon_nc
h_nc = l2grid.h_nc

# LA/OC region
region_name = 'laoc'
#lat_min = 33.5
#lat_max = 34.1
#lon_min = -118.9
#lon_max = -117.82
lat_min = 33.2
lat_max = 34.4
lon_min = -119.2
lon_max = -117.6

# full grid
#region_name = 'grid'
#lat_min = np.nanmin(lat_nc)
#lat_max = np.nanmax(lat_nc)
#lon_min = np.nanmin(lon_nc)
#lon_max = np.nanmax(lon_nc)

extent = [lon_min,lon_max,lat_min,lat_max]

figw = 18
figh = 8

axis_tick_size = 16

# large pipes lat and lon
major_nc = Dataset('/data/project1/minnaho/potw_outfall_data/updated_2013_2017/major_potw_data/major_potw_1971_2017_monthly.nc','r')
lat_potw = np.array(major_nc.variables['latitude'])
lon_potw = np.array(major_nc.variables['longitude'])

coast_10m = cpf.NaturalEarthFeature('physical','coastline','10m')

# max and min of color bar
if perc == True:
    # 1998
    v_max = 15
    v_min = -15
    # 1999
    #v_max = 50
    #v_min = -30
else:
    v_max = 5
    v_min = 0

vc_max = 1
vc_min = 0

for y in range(start_year,end_year+1):
    # if we are on the first year, starts at s_m
    if y == start_year:
        s_m = start_month
    else:
        s_m = 1
    # if we are on the last year, end at e_m
    if y == end_year:
        e_m = end_month+1
    else:
        e_m = 13
    for m in range(s_m,e_m):
        dtstr = 'Y'+str(y)+'M'+'%02d'%m
        fig,ax = plt.subplots(2,4,figsize=[figw,figh],subplot_kw=dict(projection=ccrs.PlateCarree()))
        cntrlnc = Dataset(cntrlpath+filename+dtstr+'.nc','r')
        varcnc = np.squeeze(cntrlnc.variables[var_nc])
        varcnc[varcnc>1E10] = np.nan
        varcnc[varcnc<=0] = np.nan
        varcplt = np.ones((varcnc.shape[1],varcnc.shape[2]))*np.nan
        for y_i in range(varcnc.shape[1]):
            for x_i in range(varcnc.shape[2]):
                varcplt[y_i,x_i] = np.nanmax(varcnc[:,y_i,x_i])

        p_plot1 = ax.flat[0].pcolormesh(lon_nc,lat_nc,varcplt,transform=ccrs.PlateCarree(),cmap=c_map1,vmin=vc_min,vmax=vc_max)
        ax.flat[0].set_title('CTRL',fontsize=axis_tick_size)
        # colorbar
        p0 = ax.flat[0].get_position().get_points().flatten()
        p1 = ax.flat[0].get_position().get_points().flatten()
        cb_ax1 = fig.add_axes([p0[2]-.01,p1[1]+.075,.01,p0[3]-p1[1]-.075])
        
        cb = fig.colorbar(p_plot1,cax=cb_ax1,orientation='vertical')
        #cb.set_label(cblabel,fontsize=axis_tick_size)
        cb.ax.tick_params(axis='both',which='major',labelsize=axis_tick_size-2)
        for e_i in range(len(exp)):
            print(dtstr,exp[e_i])
            datanc = Dataset(outpath+exp[e_i]+'/monthly/'+filename+dtstr+'.nc','r')
            vardnc = np.squeeze(datanc.variables[var_nc])

            vardnc[vardnc>1E10] = np.nan
            vardnc[vardnc<=0] = np.nan
            vardplt = np.ones((vardnc.shape[1],vardnc.shape[2]))*np.nan
            for y_i in range(vardnc.shape[1]):
                for x_i in range(vardnc.shape[2]):
                    vardplt[y_i,x_i] = np.nanmax(vardnc[:,y_i,x_i])

            # plot maps
            p_plot = ax.flat[e_i+1].pcolormesh(lon_nc,lat_nc,vardplt,transform=ccrs.PlateCarree(),cmap=c_map,vmin=v_min,vmax=v_max)
            
            ax.flat[e_i+1].set_title(title_exp[e_i],fontsize=axis_tick_size)
    
    
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
            if a_i >= 1 and a_i != 4:
                gl.ylabels_left = False
            if a_i >= 0 and a_i <= 3:
                gl.xlabels_bottom = False
            if region_name != 'grid':
                step_lon = .4
                step_lat = .2
            else:
                step_lon = 2
                step_lat = 1
            gl.xlocator = mticker.FixedLocator(np.arange(lon_min-step_lon,lon_max+step_lon,step_lon).astype(int))
            gl.ylocator = mticker.FixedLocator(np.arange(lat_min-step_lat,lat_max+step_lat,step_lat).astype(int))
            #gl.ylocator = mticker.FixedLocator([33.4, 33.5, 33.6, 33.7, 33.8, 33.9, 34. , 34.1, 34.2])
            #gl.ylocator = mticker.FixedLocator([33.5, 33.7, 33.9, 34.1])
            gl.yformatter = LATITUDE_FORMATTER
            gl.xformatter = LONGITUDE_FORMATTER
        
        # colorbar 
        p0 = ax.flat[3].get_position().get_points().flatten()
        p1 = ax.flat[7].get_position().get_points().flatten()
        cb_ax = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])
        
        cb = fig.colorbar(p_plot,cax=cb_ax,orientation='vertical')
        cb.set_label(cblabel,fontsize=axis_tick_size)
        cb.ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)
        
        savename = 'map_'+varstr+'_'+region_name+'_'+exp[-1]+'_'+dtstr+'.png'
        
        #plt.tight_layout()
        
        fig.savefig(savepath+savename,bbox_inches='tight')
        print(savename)
        #plt.close()

