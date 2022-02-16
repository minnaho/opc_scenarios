################################################
# plot cross section of l1617,PNDN_only,fndn90 minus CTRL 
# at major pipes
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
import matplotlib.ticker as ticker
import cmocean
import datetime as datetime
import calendar
import seawater as sw

# plot fresh vs nutrients vs control vs full
#plt.ion()

savepath = './figs/cs/'
loc = 'OCSan'

# ROMS output location
outpath = '/data/project6/minnaho/opc_scenarios/ext_depth/'

# roms var
var_name = 'O2' 
var_nc = 'var' 
cblabel = 'mmol '+var_name+' m$^{-3}$'


# scenario names 
#exp = ['PNDN_only','fndn90']
#exp = ['l1617','PNDN_only','pndn50','pndn90']
#title_exp = ['Loads 16-17','PNDN only','PNDN 50','PNDN 90']
exp = ['FNDN_only','fndn50','fndn90']
title_exp = ['FNDN only','FNDN 50','FNDN 90']

#yearmonth = 'Y1999M07_09'
#filest = 'ext_0_80_'+var_name+'_avg_'+yearmonth+'_'
#fileen = '-cntrl.nc'

yearmonth = 'summer1998'
filest = 'sub_avg_'+yearmonth+'_0_80_'+var_name+'_'
fileen = '.nc'

fpath = []
for e_i in range(len(exp)):
    fpath.append(outpath+filest+exp[e_i]+fileen)

# contour lines
if var_nc == 'rho':
    rho0 = 1027.4
    clines = [1023.5,1024,1024.5,1025,1025.5,1026,1026.5]
if var_nc == 'NH4' or var_nc == 'NO3':
    clines = [0.5,1,5,10,13]
if var_nc == 'salt':
    clines = [32.75,33,33.25,33.5,33.75,34,34.25]
if var_nc == 'temp':
    clines = [9,10,11,12,13,14,15,16,17,18,19,20]
if var_nc == 'biomass':
    clines = [0.5,1,3,5,10,15,20,25]
if var_name == 'O2':
    clines = [-25,-20,-15,-10,-5,-1,1,5,10,15,20,25]


# grid path
grid_nc = l2grid.grid_nc
lat_nc = l2grid.lat_nc
lon_nc = l2grid.lon_nc
h_nc = l2grid.h_nc

if loc == 'HTP':
    # north pipe
    lat_site = 33.920625
    lon_site = -118.528712
    ind_st = -118.65
    ind_en = -118.47
    
    dp_st = 0
    dp_en = 80

    dpipe = 60
    
    #p_loc = 551 # location to mark pipe

if loc == 'JWPCP':
    # 65% diffuser
    #lat_site = 33.6892
    #lon_site = -118.3167

    # 35% diffuser
    lat_site = 33.700737
    lon_site = -118.341962
    ind_st = -118.43
    ind_en = -118.31
    
    dp_st = 0
    dp_en = 80

    dpipe = 50

# ocsd pipe - take cross section here
if loc == 'OCSan':
    lat_site = 33.57667
    lon_site = -118.01
    #ind_st_p = 500 # remove all offshore
    #ind_en_p = 571 # remove land
    ind_st = -118.2
    ind_en = -117.963
    
    dp_st = 0
    dp_en = 80
    
    p_loc = 551 # location to mark pipe

    dpipe = 60

if loc == 'PLWTP':
    # north pipe
    lat_site = 32.671671
    lon_site = -117.325556 
    ind_st = -117.4
    ind_en = -117.25
    
    dp_st = -110
    dp_en = 0

    dpipe = 95


# calculate i and j
min_1D = np.abs( (lat_nc - lat_site)**2 + (lon_nc - lon_site)**2)
y_site, x_site = np.unravel_index(min_1D.argmin(), min_1D.shape)

# get longitudinal slice
lon_slice = lon_nc[y_site,:]
h_slice = h_nc[y_site,:]

# reshape to match z_r
# 60 by 602 in L2
# 41 by 602 in 0-80 zsliced outputs
srho_shape = 41
lon_slice_l = list(lon_slice)*srho_shape
lon_reshape = np.array(lon_slice_l).reshape(srho_shape,lon_nc.shape[1]) 

# bounds to draw contours
ind_st_p = np.nanmin(np.unique(np.where((lon_reshape[:,:]>ind_st)&(lon_reshape<ind_en))[1]))
ind_en_p = np.nanmax(np.unique(np.where((lon_reshape[:,:]>ind_st)&(lon_reshape<ind_en))[1]))

figw = 10
figh = 14

if var_nc == 'w' or var_name == 'O2':
    c_map = cmocean.cm.balance
else:
    c_map = cmocean.cm.dense

# size of marker for pipe location
psize = 300

axis_tick_size = 14

savename = 'cs_'+loc+'_avg_'+var_name+'_'+yearmonth+'.png'
fig,ax = plt.subplots(len(exp),1,figsize=[figw,figh])

z_r = np.array(Dataset(fpath[0],'r').variables['depth'])
#z_r_l = list(np.array(Dataset(fpath[0],'r').variables['depth']))*lon_slice.shape[0]
#z_r = np.array(list(np.array(Dataset(fpath[0],'r').variables['depth']))*lon_slice.shape[0]).reshape(lon_nc.shape[1],srho_shape).transpose()

for n_i in range(len(exp)):

    # roms field from output at y_site slice
    if var_nc == 'rho':
        roms_var = np.array(Dataset(fpath[n_i],'r').variables[var_nc][0,:,y_site,:])+rho0
    if var_nc == 'N2':
        # calculate N2 brunt vaisalla buoyancy frequency 
        # salt
        varsal = np.array(Dataset(fpath[n_i],'r').variables['salt'][0,:,y_site,:])
        # temp
        vartem = np.array(Dataset(fpath[n_i],'r').variables['temp'][0,:,y_site,:])
        # pressure
        varpre = sw.pres(z_r*-1,lat_site)
        roms_var = sw.bfrq(varsal,vartem,varpre,lat_site)[0]

    if var_nc == 'DIN':
        roms_var_nh4 = np.array(Dataset(fpath[n_i],'r').variables['NH4'][0,:,y_site,:])
        roms_var_nh4[roms_var_nh4>1E10] = 0
        roms_var_no3 = np.array(Dataset(fpath[n_i],'r').variables['NO3'][0,:,y_site,:])
        roms_var_no3[roms_var_no3>1E10] = 0
        roms_var_no2 = np.array(Dataset(fpath[n_i],'r').variables['NO2'][0,:,y_site,:])
        roms_var_no2[roms_var_no3>1E10] = 0

        roms_var = roms_var_nh4+roms_var_no3+roms_var_no2

    if var_nc == 'biomass':
        roms_var_dtc = np.array(Dataset(fpath[n_i],'r').variables['DIATC'][0,:,y_site,:])
        roms_var_dtc[roms_var_dtc>1E10] = 0
        roms_var_spc = np.array(Dataset(fpath[n_i],'r').variables['SPC'][0,:,y_site,:])
        roms_var_spc[roms_var_spc>1E10] = 0
        roms_var_dzc = np.array(Dataset(fpath[n_i],'r').variables['DIAZC'][0,:,y_site,:])
        roms_var_dzc[roms_var_dzc>1E10] = 0

        roms_var = roms_var_dtc+roms_var_spc+roms_var_dzc
    else:
        roms_var = np.array(Dataset(fpath[n_i],'r').variables[var_nc][0,:,y_site,:])
        roms_var[roms_var>1E10] = np.nan

    # max and min of color bar
    if var_nc == 'rho':
        v_max = 1027.5
        v_min = 1023.5
        #v_max = np.nanmax(roms_var_cntrl)
        #v_min = np.nanmin(roms_var_fulll)
    if var_nc == 'biomass':
        v_max = np.nanmax(roms_var)
        v_min = np.nanmin(roms_var)
    if var_nc == 'temp':
        #v_max = np.nanmax(roms_var_fulll)
        #v_min = np.nanmin(roms_var_cntrl)
        v_max = 20
        v_min = 10
    if var_nc == 'salt':
        v_max = 34
        v_min = 32.75
    if var_nc == 'NO3':
        v_max = 15
        #v_min = np.nanmin(roms_var_cntrl)
        v_min = 0
    if var_nc == 'NH4':
        v_max = 15
        #v_min = np.nanmin(roms_var_cntrl)
        v_min = 0
    if var_nc == 'w':
        v_max = np.nanmax(roms_var)
        v_min = -np.nanmax(roms_var)
    if var_name == 'O2':
        v_max = 25
        v_min = -25

    #  plot 
    p_plot = ax.flat[n_i].pcolor(lon_slice,z_r,roms_var,cmap=c_map,vmin=v_min,vmax=v_max)
    #p_plot = ax.flat[n_i].pcolor(lon_reshape[:,ind_st_p:ind_en_p],z_r[:,ind_st_p:ind_en_p],roms_var[:,ind_st_p:ind_en_p],cmap=c_map,vmin=v_min,vmax=v_max)

    # plot (no vmin/vmax) 
    #p_plot = ax.flat[n_i].pcolor(lon_slice,z_r,roms_var,cmap=c_map)

    # mark pipe location
    if loc == 'OCSan' or loc == 'JWPCP':
        ax.flat[n_i].scatter(lon_site+.006,dpipe,marker='^',facecolors='orange',s=psize)
    if loc == 'HTP' or loc == 'PLWTP':
        ax.flat[n_i].scatter(lon_site,dpipe,marker='^',facecolors='orange',s=psize)
    ax.flat[n_i].set_xlim([ind_st,ind_en])
    ax.flat[n_i].set_ylim([dp_st,dp_en])
    ax.flat[n_i].set_ylabel('Depth (m)',fontsize=axis_tick_size)

    # tick spacing 
    tick_spacingx = 0.05
    ax.flat[n_i].xaxis.set_major_locator(ticker.MultipleLocator(tick_spacingx))
    
    tick_spacingy = 20
    ax.flat[n_i].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacingy))

    ax.flat[n_i].set_title(title_exp[n_i],fontsize=axis_tick_size)
    ax.flat[n_i].tick_params(axis='both',which='major',labelsize=axis_tick_size)
    ax.flat[n_i].invert_yaxis()
    # tick marks on all axes
    ax.flat[n_i].yaxis.set_ticks_position('both')
    ax.flat[n_i].xaxis.set_ticks_position('both')

    # plot contours
    clinecolor = 'k'
    c_plt = ax.flat[n_i].contour(lon_reshape[0,ind_st_p:ind_en_p],z_r,roms_var[:,ind_st_p:ind_en_p],clines,colors=clinecolor,linewidths=1)

    ax.flat[n_i].clabel(c_plt,fontsize=9,fmt='%.1f',inline=1)

# outside axes loop 
p0 = ax.flat[0].get_position().get_points().flatten()
p1 = ax.flat[n_i].get_position().get_points().flatten()
cb_ax = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])

if var_nc == 'w':
    cb = fig.colorbar(p_plot,cax=cb_ax,orientation='vertical',ticks=np.linspace(v_min,v_max,6))
elif var_nc == 'O2':
    cb = fig.colorbar(p_plot,cax=cb_ax,orientation='vertical',format='%d',ticks=np.arange(v_min,v_max+25,25))
else:
    #cb = fig.colorbar(p_plot,cax=cb_ax,orientation='vertical',format='%.1f',ticks=np.linspace(v_min,v_max,11))
    cb = fig.colorbar(p_plot,cax=cb_ax,orientation='vertical',format='%.1i',ticks=np.linspace(v_min,v_max,11))

ax.flat[n_i].set_xlabel('Longitude',fontsize=axis_tick_size)
cb.set_label(cblabel,fontsize=axis_tick_size)
cb.ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)


# remove labels
for a_i in range(len(exp)-1):
    ax.flat[a_i].get_xaxis().set_ticklabels([])

#fig.suptitle('Average O2 '+loc,fontsize=axis_tick_size)


fig.savefig(savepath+savename,bbox_inches='tight')
print(savename)
#plt.close()

