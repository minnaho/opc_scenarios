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

# roms var
var_nc1 = 'NH4'
varstr1 = 'NH4'

var_nc2 = 'var'
varstr2 = 'biomass'

var_nc3 = 'O2'
varstr3 = 'O2'

# output paths 
outpath1 = '/data/project6/ROMS/L2SCB_OPC/'
outpath2 = '/data/project6/kesf/ROMS/L2SCB_AP/extract_2d/'
outpath3 = '/data/project6/ROMS/L2SCB_OPC/'

savepath = './figs/maps/'

timeunit = 'days since 1997-08-01'

#exp = ['PNDN_only','pndn50','pndn90']
#title_exp = ['PNDN only','PNDN 50','PNDN 90']
#exp = ['l1617','FNDN_only','fndn50','fndn90']
#title_exp = ['Loads 16-17','FNDN only','FNDN 50','FNDN 90']
#exp = ['l1617']
#title_exp = ['Loads 16-17']
exp = ['2013_2017']
exp2 = ['l1617']
title_exp = ['Loads 16-17']

filename1 = 'l2_scb_avg.'
filename2d = 'intbiomass_L2_2012_2017.nc'
filename2c = 'intbiomass_L2_2012_2017_nat1.nc'


# choose years
start_year = 2012
end_year = 2017

# choose months between 1 and 12
start_month = 8
end_month = 11


cblabel1 = 'mmol N m$^{-3}$'
if perc == True:
    cblabel1 = varstr2+' % change'
cblabel2 = 'mg C m$^{-3}$'
if perc == True:
    cblabel2 = varstr2+' % change'
cblabel3 = 'mmol O m$^{-3}$'
if perc == True:
    cblabel3 = varstr3+' % change'

s2d = 86400

# color map
c_map1 = cmocean.cm.dense
#c_map = cmocean.cm.thermal
c_map2 = cmocean.cm.algae
#c_map = cmocean.cm.deep
#c_map = cmocean.cm.delta
#c_map = 'PRGn'
c_map3 = cmocean.cm.ice

# control scenario
# to subtract oxygen
cntrlpath = '/data/project3/minnaho/postprocessing/cntrl_daily/'

# grid path
grid_nc = l2grid.grid_nc
lat_nc = l2grid.lat_nc
lon_nc = l2grid.lon_nc
h_nc = l2grid.h_nc

# LA/OC region
#region_name = 'laoc'
#lat_min = 33
#lat_max = 34.2
#lon_min = -119.2
#lon_max = -117.6

region_name = 'wider'
lat_min = 32
lat_max = 34.7
lon_min = -120.2
lon_max = -117

# full grid
#region_name = 'grid'
#lat_min = np.nanmin(lat_nc)
#lat_max = np.nanmax(lat_nc)
#lon_min = np.nanmin(lon_nc)
#lon_max = np.nanmax(lon_nc)

extent = [lon_min,lon_max,lat_min,lat_max]

figw = 10
figh = 6

axis_tick_size = 16

# large pipes lat and lon
major_nc = Dataset('/data/project1/minnaho/potw_outfall_data/updated_2013_2017/major_potw_data/major_potw_1971_2017_monthly.nc','r')
lat_potw = np.array(major_nc.variables['latitude'])
lon_potw = np.array(major_nc.variables['longitude'])

coast_10m = cpf.NaturalEarthFeature('physical','coastline','10m')

v_min1 = 0
v_max1 = 10

v_min2 = 1
v_max2 = 300

v_min3 = -50
v_max3 = 0

months_w_31_days = [1,3,5,7,8,10,12]
leap_years = [1992,1996,2000,2004,2008,2012,2016,2020]

varrnc2 = np.squeeze(Dataset(outpath2+filename2d,'r').variables[var_nc2])
varrnc2[varrnc2>1E10] = np.nan
varrnc2[varrnc2<=0] = np.nan

varcnc2 = np.squeeze(Dataset(outpath2+filename2c,'r').variables[var_nc2])
varcnc2[varcnc2>1E10] = np.nan
varcnc2[varcnc2<=0] = np.nan

vardplt2 = varrnc2 - varcnc2

nc_i = 0
for y_i in range(start_year,end_year+1):
    # if we are on the first year, starts at s_m
    if y_i == start_year:
        s_m = start_month
    else:
        s_m = 1
    # if we are on the last year, end at e_m
    if y_i == end_year:
        e_m = end_month+1
    else:
        e_m = 13
    for m_i in range(s_m,e_m):
        if m_i in months_w_31_days:
            ndays = 31
        if m_i not in months_w_31_days:
            ndays = 30
            if m_i == 2 and y_i in leap_years:
                ndays = 29
            if m_i == 2 and y_i not in leap_years:
                ndays = 28
        for d_i in list(range(1,ndays+1)):
            dtstr = 'Y'+str(y_i)+'M'+'%02d'%m_i+'D'+'%02d'%d_i
            print(dtstr)
            fig,ax = plt.subplots(1,1,figsize=[figw,figh],subplot_kw=dict(projection=ccrs.PlateCarree()))
            fig.suptitle(str(y_i)+'-'+'%02d'%m_i+'-'+'%02d'%d_i,fontsize=axis_tick_size)


            # plot maps
            p_plot2 = ax.pcolormesh(lon_nc,lat_nc,vardplt2[nc_i],transform=ccrs.PlateCarree(),cmap=c_map2,vmin=v_min2,vmax=v_max2,norm=mcolors.LogNorm())

            nc_i += 1
            
            ax.set_title('integrated '+varstr2,fontsize=axis_tick_size)

            # colorbar 
            p0 = ax.get_position().get_points().flatten()
            p1 = ax.get_position().get_points().flatten()
            cb_ax2 = fig.add_axes([p0[2]+.005,p1[1],.01,p0[3]-p1[1]])
            
            cb1 = fig.colorbar(p_plot2,cax=cb_ax2,orientation='vertical')
            cb1.set_label(cblabel2,fontsize=axis_tick_size-2)
            cb1.ax.tick_params(axis='both',which='major',labelsize=axis_tick_size-2)
            
            # mark pipe location
            for l_i in range(len(lon_potw)):
                ax.scatter(lon_potw[l_i],lat_potw[l_i],marker='o',facecolors='none',edgecolors='blue',s=100)
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
            if region_name == 'laoc':
                step_lon = .4
                step_lat = .3
            if region_name == 'wider':
                step_lon = 1
                step_lat = .5
            else:
                step_lon = 2
                step_lat = 1
            if region_name == 'grid': 
                gl.xlocator = mticker.FixedLocator(np.arange(lon_min-step_lon,lon_max+step_lon,step_lon).astype(int))
                gl.ylocator = mticker.FixedLocator(np.arange(lat_min-step_lat,lat_max+step_lat,step_lat).astype(int))
            else:
                gl.xlocator = mticker.FixedLocator(np.arange(lon_min-step_lon,lon_max+step_lon,step_lon).astype(int))
                gl.ylocator = mticker.FixedLocator(np.arange(lat_min-step_lat,lat_max+step_lat,step_lat))
                #gl.ylocator = mticker.FixedLocator([33.5, 33.7, 33.9, 34.1])
            gl.yformatter = LATITUDE_FORMATTER
            gl.xformatter = LONGITUDE_FORMATTER

            #fig.subplots_adjust(wspace=0.35)
            #exit()

            
            savename = 'map_vertmax_subbio_'+region_name+'_'+exp[-1]+'_'+dtstr+'.png'
            
            fig.savefig(savepath+savename,bbox_inches='tight')
            print(savename)
            plt.close()

