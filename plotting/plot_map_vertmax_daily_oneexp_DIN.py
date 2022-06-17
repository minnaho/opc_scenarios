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
var_nc1a = 'NH4'
var_nc1b = 'NO3'
varstr1 = 'DIN'

var_nc2 = 'var'
varstr2 = 'biomass'

var_nc3 = 'O2'
varstr3 = 'O2'

# output paths 
outpath1 = '/data/project6/ROMS/L2SCB_OPC/'
outpath2 = '/data/project6/minnaho/opc_scenarios/ts_int_sli/'
outpath3 = '/data/project6/ROMS/L2SCB_OPC/'

savepath = './figs/maps/'

timeunit = 'days since 1997-08-01'

#exp = ['PNDN_only','pndn50','pndn90']
#title_exp = ['PNDN only','PNDN 50','PNDN 90']
#exp = ['l1617','FNDN_only','fndn50','fndn90']
#title_exp = ['Loads 16-17','FNDN only','FNDN 50','FNDN 90']
#exp = ['l1617']
#title_exp = ['Loads 16-17']
exp = ['loads1617']
exp2 = ['l1617']
title_exp = ['Loads 16-17']

filename1 = 'l2_scb_avg.'
filename2 = 'int_100m_'+exp2[0]+'_biomass_'


# choose years
start_year = 1997
end_year = 1999

# choose months between 1 and 12
start_month = 11
end_month = 11


cblabel1 = 'mmol N m$^{-3}$'
if perc == True:
    cblabel1 = varstr2+' % change'
cblabel2 = 'mmol C m$^{-3}$'
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
region_name = 'laoc'
#lat_min = 33.5
#lat_max = 34.1
#lon_min = -118.9
#lon_max = -117.82
lat_min = 33
lat_max = 34.2
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

v_min1 = 0
v_max1 = 10

v_min2 = 100
v_max2 = 1000

v_min3 = -50
v_max3 = 0

vc_max = 1
vc_min = 0

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
        if m_i in months_w_31_days:
            ndays = 31
        if m_i not in months_w_31_days:
            ndays = 30
            if m_i == 2 and y_i in leap_years:
                ndays = 29
            if m_i == 2 and y_i not in leap_years:
                ndays = 28
        if y_i == 1999 and m_i == 11:
            ndays = 27 # last day of the simulation
        for d_i in list(range(1,ndays+1)):
            dtstr = 'Y'+str(y_i)+'M'+'%02d'%m_i+'D'+'%02d'%d_i
            print(dtstr,exp[0])
            fig,ax = plt.subplots(1,3,figsize=[figw,figh],subplot_kw=dict(projection=ccrs.PlateCarree()))
            cntrlnc = Dataset(cntrlpath+filename1+dtstr+'.nc','r')
            varcnc = np.squeeze(cntrlnc.variables[var_nc3])
            varcnc[varcnc>1E10] = np.nan
            varcnc[varcnc<=0] = np.nan

            # var1 - NH4
            datanc = Dataset(outpath1+exp[0]+'/daily/'+filename1+dtstr+'.nc','r')
            vardnc1a = np.squeeze(datanc.variables[var_nc1a])
            vardnc1b = np.squeeze(datanc.variables[var_nc1b])

            vardnc1a[vardnc1a>1E10] = np.nan
            vardnc1a[vardnc1a<=0] = np.nan
            vardnc1b[vardnc1b>1E10] = np.nan
            vardnc1b[vardnc1b<=0] = np.nan

            vardnc1 = vardnc1a+vardnc1b

            # get vertical maximum
            vardplt1 = np.ones((vardnc1.shape[1],vardnc1.shape[2]))*np.nan
            for eta_i in range(vardnc1.shape[1]):
                for xi_i in range(vardnc1.shape[2]):
                    vardplt1[eta_i,xi_i] = np.nanmax(vardnc1[:,eta_i,xi_i])
            # var2 - biomass (already have sum of biomass)
            vardplt2 = np.squeeze(Dataset(outpath2+filename2+dtstr+'.nc','r').variables[var_nc2])
            vardplt2[vardplt2>1E10] = np.nan
            vardplt2[vardplt2<=0] = np.nan
            
            # var3 - O2
            vardnc3rd = np.squeeze(datanc.variables[var_nc3])

            vardnc3rd[vardnc3rd>1E10] = np.nan
            vardnc3rd[vardnc3rd<=0] = np.nan
            vardnc3 = vardnc3rd - varcnc 
            # get vertical minimum
            vardplt3 = np.ones((vardnc3.shape[1],vardnc3.shape[2]))*np.nan
            for eta_i in range(vardnc3.shape[1]):
                for xi_i in range(vardnc3.shape[2]):
                    vardplt3[eta_i,xi_i] = np.nanmin(vardnc3[:,eta_i,xi_i])

            # plot maps
            p_plot1 = ax.flat[0].pcolormesh(lon_nc,lat_nc,vardplt1,transform=ccrs.PlateCarree(),cmap=c_map1,vmin=v_min1,vmax=v_max1)
            p_plot2 = ax.flat[1].pcolormesh(lon_nc,lat_nc,vardplt2,transform=ccrs.PlateCarree(),cmap=c_map2,vmin=v_min2,vmax=v_max2,norm=mcolors.LogNorm())
            p_plot3 = ax.flat[2].pcolormesh(lon_nc,lat_nc,vardplt3,transform=ccrs.PlateCarree(),cmap=c_map3,vmin=v_min3,vmax=v_max3)
            
            ax.flat[0].set_title(varstr1+' maximum',fontsize=axis_tick_size)
            ax.flat[1].set_title('integrated '+varstr2,fontsize=axis_tick_size)
            ax.flat[2].set_title(varstr3+' minimum',fontsize=axis_tick_size)

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
            # colorbar 
            p0 = ax.flat[2].get_position().get_points().flatten()
            p1 = ax.flat[2].get_position().get_points().flatten()
            cb_ax3 = fig.add_axes([p0[2]+.005,p1[1],.01,p0[3]-p1[1]])
            
            cb3 = fig.colorbar(p_plot3,cax=cb_ax3,orientation='vertical')
            cb3.set_label(cblabel3,fontsize=axis_tick_size-2)
            cb3.ax.tick_params(axis='both',which='major',labelsize=axis_tick_size-2)
            
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
                if region_name != 'grid':
                    step_lon = .4
                    step_lat = .3
                else:
                    step_lon = 2
                    step_lat = 1
                if region_name == 'grid': 
                    gl.xlocator = mticker.FixedLocator(np.arange(lon_min-step_lon,lon_max+step_lon,step_lon).astype(int))
                    gl.ylocator = mticker.FixedLocator(np.arange(lat_min-step_lat,lat_max+step_lat,step_lat).astype(int))
                else:
                    gl.xlocator = mticker.FixedLocator(np.arange(lon_min-step_lon,lon_max+step_lon,step_lon))
                    gl.ylocator = mticker.FixedLocator(np.arange(lat_min-step_lat,lat_max+step_lat,step_lat))
                    #gl.ylocator = mticker.FixedLocator([33.5, 33.7, 33.9, 34.1])
                gl.yformatter = LATITUDE_FORMATTER
                gl.xformatter = LONGITUDE_FORMATTER

            fig.subplots_adjust(wspace=0.35)
            
            savename = 'map_vertmax_DIN_'+region_name+'_'+exp[-1]+'_'+dtstr+'.png'
            
            fig.savefig(savepath+savename,bbox_inches='tight')
            print(savename)
            plt.close()

