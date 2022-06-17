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
var_nc1 = 'var'
varstr1 = 'biomass'

var_nc2 = 'O2'
varstr2 = 'O2'

# output paths 
outpath1 = '/data/project6/minnaho/L2SCB_2012_2017/ts_int_sli/'
outpath2 = '/data/project6/ROMS/L2SCB_AP/monthly_2012_2017/'

savepath = './figs/maps/'

timeunit = 'days since 1997-08-01'

exp = ['2013_2017']

# r = read/anth, c = control
filename1r = 'int_100m_anth_biomass_'
filename1c = 'int_100m_cntrl_biomass_'
filename2 = 'l2_scb_avg.'


# choose years
start_year = 2015
end_year = 2017

# choose months between 1 and 12
start_month = 6
end_month = 11


cblabel1 = 'mmol C m$^{-3}$'
if perc == True:
    cblabel1 = varstr2+' % change'
cblabel2 = 'mmol O m$^{-3}$'
if perc == True:
    cblabel2 = varstr2+' % change'

# color map
c_map1 = cmocean.cm.algae
c_map2 = cmocean.cm.ice

# control scenario
# to subtract oxygen
cntrlpath1 = '/data/project6/minnaho/L2SCB_2012_2017/ts_int_sli/'
cntrlpath2 = '/data/project6/ROMS/L2SCB/monthly_2012_2017/'

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
lat_min = 33
lat_max = 34.2
lon_min = -119.2
lon_max = -117.6

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

figw = 18
figh = 6

axis_tick_size = 16

# large pipes lat and lon
major_nc = Dataset('/data/project1/minnaho/potw_outfall_data/updated_2013_2017/major_potw_data/major_potw_1971_2017_monthly.nc','r')
lat_potw = np.array(major_nc.variables['latitude'])
lon_potw = np.array(major_nc.variables['longitude'])

coast_10m = cpf.NaturalEarthFeature('physical','coastline','10m')

v_min1 = 100
v_max1 = 1000

v_min2 = -30
v_max2 = 0

months_w_31_days = [1,3,5,7,8,10,12]
leap_years = [1992,1996,2000,2004,2008,2012,2016,2020]

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
        dtstr = 'Y'+str(y_i)+'M'+'%02d'%m_i
        print(dtstr)

        fig,ax = plt.subplots(1,2,figsize=[figw,figh],subplot_kw=dict(projection=ccrs.PlateCarree()))
        fig.suptitle(str(y_i)+'-'+'%02d'%m_i,fontsize=axis_tick_size)

        cntrlnc1 = Dataset(cntrlpath1+filename1c+dtstr+'.nc','r')
        cntrlnc2 = Dataset(cntrlpath2+filename2+dtstr+'.nc','r')

        varcnc1 = np.squeeze(cntrlnc1.variables[var_nc1])
        varcnc1[varcnc1>1E10] = np.nan
        varcnc1[varcnc1<=0] = np.nan

        varcnc2 = np.squeeze(cntrlnc2.variables[var_nc2])
        varcnc2[varcnc2>2E20] = np.nan
        varcnc2[varcnc2<=0] = np.nan

        # var1 - biomass
        varrnc1 = np.squeeze(Dataset(outpath1+filename1r+dtstr+'.nc','r').variables[var_nc1])
        varrnc1[varrnc1>1E10] = np.nan
        varrnc1[varrnc1<=0] = np.nan

        vardplt1 = varrnc1

        # var2 - O2
        vardnc2rd = np.squeeze(Dataset(outpath2+filename2+dtstr+'.nc','r').variables[var_nc2])

        vardnc2rd[vardnc2rd>1E10] = np.nan
        vardnc2rd[vardnc2rd<=0] = np.nan
        vardnc2 = vardnc2rd - varcnc2
        # get vertical minimum
        #vardplt2 = np.nanmin(vardnc2,axis=0)
        # vertical lowest 10th percentile
        vardplt2 = np.nanpercentile(vardnc2,10,axis=0)

        # plot maps
        p_plot1 = ax.flat[0].pcolormesh(lon_nc,lat_nc,vardplt1,transform=ccrs.PlateCarree(),cmap=c_map1,vmin=v_min1,vmax=v_max1,norm=mcolors.LogNorm())
        p_plot2 = ax.flat[1].pcolormesh(lon_nc,lat_nc,vardplt2,transform=ccrs.PlateCarree(),cmap=c_map2,vmin=v_min2,vmax=v_max2)
        
        ax.flat[0].set_title('integrated '+varstr1,fontsize=axis_tick_size)
        ax.flat[1].set_title(varstr2+' minimum',fontsize=axis_tick_size)

        # colorbar 
        p0 = ax.flat[0].get_position().get_points().flatten()
        p1 = ax.flat[0].get_position().get_points().flatten()
        cb_ax1 = fig.add_axes([p0[2]+.005,p1[1],.01,p0[3]-p1[1]])
        
        cb1 = fig.colorbar(p_plot1,cax=cb_ax1,orientation='vertical')
        cb1.set_label(cblabel1,fontsize=axis_tick_size-2)
        cb1.ax.tick_params(axis='both',which='major',labelsize=axis_tick_size-2)
        # colorbar 
        p0 = ax.flat[1].get_position().get_points().flatten()
        p1 = ax.flat[1].get_position().get_points().flatten()
        cb_ax2 = fig.add_axes([p0[2]+.005,p1[1],.01,p0[3]-p1[1]])
        
        cb2 = fig.colorbar(p_plot2,cax=cb_ax2,orientation='vertical')
        cb2.set_label(cblabel2,fontsize=axis_tick_size-2)
        cb2.ax.tick_params(axis='both',which='major',labelsize=axis_tick_size-2)

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
        
        savename = 'map_vertmax_monthly_totbio_O2perc10_'+region_name+'_'+exp[-1]+'_'+dtstr+'.png'
        
        fig.savefig(savepath+savename,bbox_inches='tight')
        print(savename)
        plt.close()

