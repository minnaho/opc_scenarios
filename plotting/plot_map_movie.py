################################################
# plot daily maps of freshwater/nutrient/control/full
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


# plot fresh vs nutrients vs control vs full
#plt.ion()

# map, int, or sli
filename = 'sli'

savepath = './figs/map/'
depth = '40m'

timeunit = 'days since 1999-07-01'

# roms var
var_nc = 'salt'
#cblabel = 'Density (kg m$^{-3}$)'
#cblabel = 'Temperature C'
cblabel = 'PSU'
#cblabel = var_nc+' (mmol m$^{-3}$)'
#cblabel = var_nc+' (PSU)'
#cblabel = 'integrated mmol C m$^{-2}$'

# color map
c_map = cmocean.cm.dense
#c_map = cmocean.cm.thermal
#c_map = cmocean.cm.algae
#c_map = cmocean.cm.deep

# path of outputs
# 40 m slice z slice
freshnc = filename+'_'+depth+'_fresh_all.nc'
nutrinc = filename+'_'+depth+'_nutri_all.nc'
cntrlnc = filename+'_'+depth+'_cntrl_all.nc'
fulllnc = filename+'_'+depth+'_pipes_all.nc' # full is actually pipes

# other extractions from matlab
#freshnc = filename+'_'+depth+'_fresh_biomass.nc'
#nutrinc = filename+'_'+depth+'_nutri_biomass.nc'
#cntrlnc = filename+'_'+depth+'_cntrl_biomass.nc'
#fulllnc = filename+'_'+depth+'_pipes_biomass.nc' # full is actually pipes

# 40 m slice z slice
freshpath = '/data/project3/minnaho/freshwater/postprocessing/depth'+depth+'_sli/fresh/'
nutripath = '/data/project3/minnaho/freshwater/postprocessing/depth'+depth+'_sli/nutri/'
cntrlpath = '/data/project3/minnaho/freshwater/postprocessing/depth'+depth+'_sli/cntrl/'
fulllpath = '/data/project3/minnaho/freshwater/postprocessing/depth'+depth+'_sli/pipes/'

# other extractions from matlab
#freshpath = '/data/project3/minnaho/freshwater/postprocessing/ts_int_sli/'
#nutripath = '/data/project3/minnaho/freshwater/postprocessing/ts_int_sli/'
#cntrlpath = '/data/project3/minnaho/freshwater/postprocessing/ts_int_sli/'
#fulllpath = '/data/project3/minnaho/freshwater/postprocessing/ts_int_sli/'

freshnc = Dataset(freshpath+freshnc,'r')
nutrinc = Dataset(nutripath+nutrinc,'r')
cntrlnc = Dataset(cntrlpath+cntrlnc,'r')
fulllnc = Dataset(fulllpath+fulllnc,'r')

dtshape = freshnc.variables[var_nc].shape[0]
dt = num2date(range(dtshape),timeunit,only_use_cftime_datetimes=False,only_use_python_datetimes=True)

# grid path
grid_nc = l2grid.grid_nc
lat_nc = l2grid.lat_nc
lon_nc = l2grid.lon_nc
h_nc = l2grid.h_nc

lat_min = 33.5
lat_max = 34.1
lon_min = -118.9
lon_max = -117.82

extent = [lon_min,lon_max,lat_min,lat_max]

figw = 12
figh = 7

axis_tick_size = 16

#HTP
# north pipe
lat_site_htp = 33.920625 # pipe coordinates as implemented in model
lon_site_htp = -118.529712
#JWPCP
# 35% diffuser
lat_site_jwp = 33.700737
lon_site_jwp = -118.341962
#OCSD
lat_site_ocs = 33.57667
lon_site_ocs = -118.01
#PLWTP
# north pipe
lat_site_plw = 32.671671
lon_site_plw = -117.325556

lat_potw = [lat_site_htp,lat_site_jwp,lat_site_ocs,lat_site_plw]
lon_potw = [lon_site_htp,lon_site_jwp,lon_site_ocs,lon_site_plw]

coast_10m = cpf.NaturalEarthFeature('physical','coastline','10m')
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
    v_max = 15
    #v_min = np.nanmin(roms_var_cntrl)
    v_min = 0
if var_nc == 'temp':
    v_max = 20
    #v_min = np.nanmin(roms_var_cntrl)
    v_min = 10
if var_nc == 'var':
    v_max = 1000
    #v_min = np.nanmin(roms_var_cntrl)
    v_min = 1


for t_i in range(len(dt)):
    fig,ax = plt.subplots(2,2,figsize=[figw,figh],subplot_kw=dict(projection=ccrs.PlateCarree()))
    roms_var_cntrl = np.array(cntrlnc.variables[var_nc][t_i,:,:])
    roms_var_fresh = np.array(freshnc.variables[var_nc][t_i,:,:])
    roms_var_nutri = np.array(nutrinc.variables[var_nc][t_i,:,:])
    roms_var_fulll = np.array(fulllnc.variables[var_nc][t_i,:,:])
    
    roms_var_cntrl[roms_var_cntrl>1E10] = np.nan
    roms_var_fresh[roms_var_fresh>1E10] = np.nan
    roms_var_nutri[roms_var_nutri>1E10] = np.nan
    roms_var_fulll[roms_var_fulll>1E10] = np.nan
    
    roms_var_cntrl[roms_var_cntrl==0] = np.nan
    roms_var_fresh[roms_var_fresh==0] = np.nan
    roms_var_nutri[roms_var_nutri==0] = np.nan
    roms_var_fulll[roms_var_fulll==0] = np.nan
    
    
    # plot maps
    #p_plot_cntrl = ax.flat[0].pcolormesh(lon_nc,lat_nc,roms_var_cntrl,transform=ccrs.PlateCarree(),cmap=c_map,vmin=v_min,vmax=v_max,norm=mcolors.LogNorm())
    #p_plot_fresh = ax.flat[1].pcolormesh(lon_nc,lat_nc,roms_var_fresh,transform=ccrs.PlateCarree(),cmap=c_map,vmin=v_min,vmax=v_max,norm=mcolors.LogNorm())
    #p_plot_nutri = ax.flat[2].pcolormesh(lon_nc,lat_nc,roms_var_nutri,transform=ccrs.PlateCarree(),cmap=c_map,vmin=v_min,vmax=v_max,norm=mcolors.LogNorm())
    #p_plot_fulll = ax.flat[3].pcolormesh(lon_nc,lat_nc,roms_var_fulll,transform=ccrs.PlateCarree(),cmap=c_map,vmin=v_min,vmax=v_max,norm=mcolors.LogNorm())
    p_plot_cntrl = ax.flat[0].pcolormesh(lon_nc,lat_nc,roms_var_cntrl,transform=ccrs.PlateCarree(),cmap=c_map,vmin=v_min,vmax=v_max)
    p_plot_fresh = ax.flat[1].pcolormesh(lon_nc,lat_nc,roms_var_fresh,transform=ccrs.PlateCarree(),cmap=c_map,vmin=v_min,vmax=v_max)
    p_plot_nutri = ax.flat[2].pcolormesh(lon_nc,lat_nc,roms_var_nutri,transform=ccrs.PlateCarree(),cmap=c_map,vmin=v_min,vmax=v_max)
    p_plot_fulll = ax.flat[3].pcolormesh(lon_nc,lat_nc,roms_var_fulll,transform=ccrs.PlateCarree(),cmap=c_map,vmin=v_min,vmax=v_max)
    
    
    # plot contours
    #clinecolor = 'k'
    #c_plt_cntrl = ax.flat[0].contour(lon_reshape[:,ind_st_p:ind_en_p],z_r_cntrl[:,ind_st_p:ind_en_p],roms_var_cntrl[:,ind_st_p:ind_en_p],clines,colors=clinecolor,linewidths=1)
    #c_plt_fresh = ax.flat[1].contour(lon_reshape[:,ind_st_p:ind_en_p],z_r_fresh[:,ind_st_p:ind_en_p],roms_var_fresh[:,ind_st_p:ind_en_p],clines,colors=clinecolor,linewidths=1)
    #c_plt_nutri = ax.flat[2].contour(lon_reshape[:,ind_st_p:ind_en_p],z_r_nutri[:,ind_st_p:ind_en_p],roms_var_nutri[:,ind_st_p:ind_en_p],clines,colors=clinecolor,linewidths=1)
    #c_plt_fulll = ax.flat[3].contour(lon_reshape[:,ind_st_p:ind_en_p],z_r_fulll[:,ind_st_p:ind_en_p],roms_var_fulll[:,ind_st_p:ind_en_p],clines,colors=clinecolor,linewidths=1)
    #
    #ax.flat[0].clabel(c_plt_cntrl,fontsize=9,fmt='%.1f',inline=1)
    #ax.flat[1].clabel(c_plt_fresh,fontsize=9,fmt='%.1f',inline=1)
    #ax.flat[2].clabel(c_plt_nutri,fontsize=9,fmt='%.1f',inline=1)
    #ax.flat[3].clabel(c_plt_fulll,fontsize=9,fmt='%.1f',inline=1)
    
    fig.suptitle(depth+' '+var_nc+' '+dt[t_i].strftime('%Y-%m-%d'),fontsize=axis_tick_size)
    
    ax.flat[0].set_title('CTRL',fontsize=axis_tick_size)
    ax.flat[1].set_title('FRESH',fontsize=axis_tick_size)
    ax.flat[2].set_title('NUTR',fontsize=axis_tick_size)
    ax.flat[3].set_title('PIPES',fontsize=axis_tick_size)
    
    
    for a_i in range(len(ax.flat)):
        # mark pipe location
        for l_i in range(len(lon_potw)):
            ax.flat[a_i].scatter(lon_potw[l_i],lat_potw[l_i],marker='o',facecolors='none',edgecolors='orange',s=100)
        # other grid stuff
        ax.flat[a_i].tick_params(axis='both',which='major',labelsize=axis_tick_size)
        ax.flat[a_i].yaxis.set_ticks_position('both')
        ax.flat[a_i].xaxis.set_ticks_position('both')
        ax.flat[a_i].add_feature(coast_10m,facecolor='None',edgecolor='k')
        ax.flat[a_i].set_extent(extent)
        # lat/lon axes
        gl = ax.flat[a_i].gridlines(draw_labels=True,linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlabel_style = {'size':axis_tick_size}
        gl.ylabel_style = {'size':axis_tick_size}
        if a_i == 0:
            gl.xlabels_bottom = False
        if a_i == 1:
            gl.xlabels_bottom = False
            gl.ylabels_left = False
        if a_i == 3:
            gl.ylabels_left = False
        step_lon = .3
        step_lat = .1
        gl.xlocator = mticker.FixedLocator(np.arange(lon_min-step_lon,lon_max+step_lon,step_lon))
        #gl.ylocator = mticker.FixedLocator(np.arange(lat_min-step_lat,lat_max+step_lat,step_lat))
        gl.ylocator = mticker.FixedLocator([33.4, 33.5, 33.6, 33.7, 33.8, 33.9, 34. , 34.1, 34.2])
        gl.yformatter = LATITUDE_FORMATTER
        gl.xformatter = LONGITUDE_FORMATTER
    
    # colorbar
    p0 = ax.flat[1].get_position().get_points().flatten()
    p1 = ax.flat[3].get_position().get_points().flatten()
    cb_ax = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])
    
    cb = fig.colorbar(p_plot_cntrl,cax=cb_ax,orientation='vertical',format='%.1f')
    cb.set_label(cblabel,fontsize=axis_tick_size)
    cb.ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)
    
    savename = 'map_'+var_nc+'_'+depth+'_'+dt[t_i].strftime('Y%YM%mD%d')
    #savename = 'map_biomass_'+depth+'_'+dt[t_i].strftime('Y%YM%mD%d')

    plt.tight_layout()
    
    fig.savefig(savepath+savename,bbox_inches='tight')
    print(savename)
    plt.close()

