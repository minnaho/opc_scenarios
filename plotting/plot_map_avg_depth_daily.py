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

# avg maps
outpath = '/data/project6/minnaho/opc_scenarios/avg_depth/'
filename = 'avg'

savepath = './figs/maps/'
depth = '30_45'

# month or season

start_year = 1997
end_year = 1998

start_month = 11
end_month = 6

exp = ['cntrl','l1617','PNDN_only','FNDN_only',
       'pndn50','fndn90']

title_exp = ['CTRL','Loads 16-17','PNDN only','FNDN only',
             'PNDN 50','FNDN 90']

# roms var
var_nc = 'NH4'
#cblabel = 'Density (kg m$^{-3}$)'
#cblabel = 'Temperature C'
cblabel = var_nc+' mmol m$^{-3}$' 
#cblabel = var_nc+' (PSU)'
#cblabel = 'integrated mmol C m$^{-2}$'

# color map
c_map = cmocean.cm.dense
#c_map = cmocean.cm.thermal
#c_map = cmocean.cm.algae
#c_map = cmocean.cm.deep

# grid path
grid_nc = l2grid.grid_nc
lat_nc = l2grid.lat_nc
lon_nc = l2grid.lon_nc
h_nc = l2grid.h_nc

# LA/OC region
lat_min = 33.5
lat_max = 34.1
lon_min = -118.9
lon_max = -117.82

#lat_min = np.nanmin(lat_nc)
#lat_max = np.nanmax(lat_nc)
#lon_min = np.nanmin(lon_nc)
#lon_max = np.nanmax(lon_nc)

extent = [lon_min,lon_max,lat_min,lat_max]

figw = 16
figh = 5

axis_tick_size = 16

# large pipes lat and lon
major_nc = Dataset('/data/project1/minnaho/potw_outfall_data/updated_2013_2017/major_potw_data/major_potw_1971_2017_monthly.nc','r')
lat_potw = np.array(major_nc.variables['latitude'])
lon_potw = np.array(major_nc.variables['longitude'])

coast_10m = cpf.NaturalEarthFeature('physical','coastline','10m')

# max and min of color bar
if var_nc == 'rho':
    v_max = 1027.5
    v_min = 1023.5
    #v_max = np.nanmax(roms_var_cntrl)
    #v_min = np.nanmin(roms_var_fulll)
if var_nc == 'salt':
    v_max = 34
    v_min = 33.1
if var_nc == 'NO3':
    v_max = 15
    #v_min = np.nanmin(roms_var_cntrl)
    v_min = 0
if var_nc == 'NH4':
    v_max = 5
    #v_min = np.nanmin(roms_var_cntrl)
    v_min = 0
if var_nc == 'temp':
    v_max = 20
    #v_min = np.nanmin(roms_var_cntrl)
    v_min = 10

months_w_31_days = [1,3,5,7,8,10,12]
leap_years = [1992,1996,2000,2004,2008,2012,2016,2020]


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
        if m in months_w_31_days:
            ndays = 31
        if m not in months_w_31_days:
            ndays = 30
            if m == 2 and y in leap_years:
                ndays = 29
            if m == 2 and y not in leap_years:
                ndays = 28
        for d in list(range(1,ndays+1)):
            year_month = 'Y'+str(y)+'M'+'%02d'%m+'D'+'%02d'%d
            fig,ax = plt.subplots(2,3,figsize=[figw,figh],subplot_kw=dict(projection=ccrs.PlateCarree()))

            for e_i in range(len(exp)):
                datanc = Dataset(outpath+filename+'_'+depth+'_'+var_nc+'_'+year_month+'_'+exp[e_i]+'.nc','r')
                varplt = np.squeeze(np.array(datanc['var']))
                varplt[varplt>1E10] = np.nan
                varplt[varplt<=0] = np.nan
    
                # plot maps
                p_plot = ax.flat[e_i].pcolormesh(lon_nc,lat_nc,varplt,transform=ccrs.PlateCarree(),cmap=c_map,vmin=v_min,vmax=v_max)
                ax.flat[e_i].set_title(title_exp[e_i],fontsize=axis_tick_size)
                # plot pipe locations
                for l_i in range(len(lon_potw)):
                    ax.flat[e_i].scatter(lon_potw[l_i],lat_potw[l_i],marker='o',facecolors='none',edgecolors='orange',s=100)

                # other grid stuff
                ax.flat[e_i].tick_params(axis='both',which='major',labelsize=axis_tick_size)
                ax.flat[e_i].yaxis.set_ticks_position('both')
                ax.flat[e_i].xaxis.set_ticks_position('both')
                ax.flat[e_i].add_feature(coast_10m,facecolor='None',edgecolor='k')
                ax.flat[e_i].set_extent(extent)
                # lat/lon axes
                gl = ax.flat[e_i].gridlines(draw_labels=True,linestyle='--')
                gl.xlabels_top = False
                gl.ylabels_right = False
                gl.xlabel_style = {'size':axis_tick_size}
                gl.ylabel_style = {'size':axis_tick_size}
                if e_i == 0:
                    gl.xlabels_bottom = False
                if e_i == 1 or e_i == 2:
                    gl.xlabels_bottom = False
                    gl.ylabels_left = False
                if e_i == 4 or e_i == 5 or e_i == 6:
                    gl.ylabels_left = False
                step_lon = .4
                step_lat = .2
                gl.xlocator = mticker.FixedLocator(np.arange(lon_min-step_lon,lon_max+step_lon,step_lon))
                #gl.ylocator = mticker.FixedLocator(np.arange(lat_min-step_lat,lat_max+step_lat,step_lat))
                #gl.ylocator = mticker.FixedLocator([33.4, 33.5, 33.6, 33.7, 33.8, 33.9, 34. , 34.1, 34.2])
                gl.ylocator = mticker.FixedLocator([33.5, 33.7, 33.9, 34.1])
                gl.yformatter = LATITUDE_FORMATTER
                gl.xformatter = LONGITUDE_FORMATTER

            # colorbar
            p0 = ax.flat[2].get_position().get_points().flatten()
            p1 = ax.flat[5].get_position().get_points().flatten()
            cb_ax = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])
            
            cb = fig.colorbar(p_plot,cax=cb_ax,orientation='vertical',format='%.1f')
            cb.set_label(cblabel,fontsize=axis_tick_size)
            cb.ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)
            
            fig.suptitle('Average NH4 30-45 m depth '+str(y)+'-'+'%02d'%m+'-'+'%02d'%d,fontsize=axis_tick_size)

            savename = 'map6_'+var_nc+'_'+depth+'_'+year_month
            
            plt.tight_layout()
            
            fig.savefig(savepath+savename,bbox_inches='tight')
            print(savename)
            plt.close()
    
    
    
