################################################
# plot map of freshwater/nutrient/control/full
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

# plot fresh vs nutrients vs control vs full
#plt.ion()

savepath = './figs/cs/'
loc = 'HTP'

# roms var
var_nc = 'temp'
#cblabel = 'Density (kg m$^{-3}$)'
cblabel = 'Temperature C'
#cblabel = 'NO3 (mmol m$^{-3}$)'
#cblabel = 'NH4 (mmol m$^{-3}$)'
#cblabel = 'salt (PSU)'

# path of outputs
freshpath = '/data/project6/ROMS/L2SCB_AP/fresh/monthly/'
nutripath = '/data/project6/ROMS/L2SCB_AP/nutrients/monthly/'
cntrlpath = '/data/project6/ROMS/L2SCB_1997_2000/monthly/'
fulllpath = '/data/project6/ROMS/L2SCB_AP/monthly/'

# choose year and month
#year = 1999
#month = 9

start_year = 1999
end_year = 2000

# between 1 and 12
start_month = 7
end_month = 9


if var_nc == 'rho':
    rho0 = 1027.4
    clines = [1023.5,1024,1024.5,1025,1025.5,1026,1026.5]
if var_nc == 'NH4' or var_nc == 'NO3':
    clines = [0.5,1,5,10,13]
if var_nc == 'salt':
    clines = [32.75,33,33.25,33.5,33.75,34,34.25]
if var_nc == 'temp':
    clines = [9,10,11,12,13,14,15,16,17,18,19,20]

#savename = 'cs_'+var_nc+'_Y'+str(year)+'M'+'%02d'%month+'.png'

# roms file name
#ncfile = 'l2_scb_avg.Y'+str(year)+'M'+'%02d'%month+'.nc'

ncfile = []
savename = []
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
        year_month = 'Y'+str(y)+'M'+'%02d'%m
        ncfile.append('l2_scb_avg.'+year_month+'.nc')
        savename.append('cs_'+loc+'_'+var_nc+'_Y'+str(y)+'M'+'%02d'%m+'.png')

# outputs
freshnc = []
nutrinc = []
cntrlnc = []
fulllnc = []

for n_i in range(len(ncfile)):
    freshnc.append(Dataset(freshpath+ncfile[n_i],'r'))
    nutrinc.append(Dataset(nutripath+ncfile[n_i],'r'))
    cntrlnc.append(Dataset(cntrlpath+ncfile[n_i],'r'))
    fulllnc.append(Dataset(fulllpath+ncfile[n_i],'r'))


# grid path
grid_nc = l2grid.grid_nc
lat_nc = l2grid.lat_nc
lon_nc = l2grid.lon_nc
h_nc = l2grid.h_nc

if loc == 'HTP':
    # north pipe
    lat_site = 33.920625
    lon_site = -118.529712
    ind_st = -118.6
    ind_en = -118.46
    
    dp_st = -80
    dp_en = 0

    dpipe = -60
    
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
    
    dp_st = -80
    dp_en = 0

    dpipe = -50

# ocsd pipe - take cross section here
if loc == 'OCSD':
    lat_site = 33.57667
    lon_site = -118.01
    #ind_st_p = 500 # remove all offshore
    #ind_en_p = 571 # remove land
    ind_st = -118.1
    ind_en = -117.963
    
    dp_st = -70
    dp_en = 0
    
    p_loc = 551 # location to mark pipe

if loc == 'PLWTP':
    # north pipe
    lat_site = 32.671671
    lon_site = -117.325556 
    ind_st = -117.4
    ind_en = -117.25
    
    dp_st = -110
    dp_en = 0

    dpipe = -95


# calculate i and j
min_1D = np.abs( (lat_nc - lat_site)**2 + (lon_nc - lon_site)**2)
y_site, x_site = np.unravel_index(min_1D.argmin(), min_1D.shape)

# get longitudinal slice
lon_slice = lon_nc[y_site,:]
h_slice = h_nc[y_site,:]

# reshape to match z_r
# 60 by 602 in L2
lon_slice_l = list(lon_slice)*freshnc[0].variables['temp'].shape[1]
lon_reshape = np.array(lon_slice_l).reshape(freshnc[0].variables['temp'].shape[1],lon_nc.shape[1]) 

# bounds to draw contours
ind_st_p = np.nanmin(np.unique(np.where((lon_reshape[:,:]>ind_st)&(lon_reshape<ind_en))[1]))
ind_en_p = np.nanmax(np.unique(np.where((lon_reshape[:,:]>ind_st)&(lon_reshape<ind_en))[1]))

figw = 14
figh = 7.5
c_map = cmocean.cm.dense

msize = 50
psize = 300

axis_tick_size = 14

for n_i in range(len(ncfile)):
    fig,ax = plt.subplots(2,2,figsize=[figw,figh])
    # roms field from output at y_site slice
    if var_nc == 'rho':
        roms_var_cntrl = np.array(cntrlnc[n_i].variables[var_nc][0,:,y_site,:])+rho0
        roms_var_fresh = np.array(freshnc[n_i].variables[var_nc][0,:,y_site,:])+rho0
        roms_var_nutri = np.array(nutrinc[n_i].variables[var_nc][0,:,y_site,:])+rho0
        roms_var_fulll = np.array(fulllnc[n_i].variables[var_nc][0,:,y_site,:])+rho0
    else:
        roms_var_cntrl = np.array(cntrlnc[n_i].variables[var_nc][0,:,y_site,:])
        roms_var_fresh = np.array(freshnc[n_i].variables[var_nc][0,:,y_site,:])
        roms_var_nutri = np.array(nutrinc[n_i].variables[var_nc][0,:,y_site,:])
        roms_var_fulll = np.array(fulllnc[n_i].variables[var_nc][0,:,y_site,:])
    
    roms_var_cntrl[roms_var_cntrl>1E10] = np.nan
    roms_var_fresh[roms_var_fresh>1E10] = np.nan
    roms_var_nutri[roms_var_nutri>1E10] = np.nan
    roms_var_fulll[roms_var_fulll>1E10] = np.nan
    
    # max and min of color bar
    if var_nc == 'rho':
        v_max = 1027.5
        v_min = 1023.5
        #v_max = np.nanmax(roms_var_cntrl)
        #v_min = np.nanmin(roms_var_fulll)
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
    
    
    # get depths at this y_site slice for each scenario
    z_r_cntrl = depths.get_zr_zw_tind(cntrlnc[n_i],grid_nc,0,[y_site-1,y_site+1,0,cntrlnc[n_i].variables[var_nc].shape[3]])[0][:,1,:]
    z_r_fresh = depths.get_zr_zw_tind(freshnc[n_i],grid_nc,0,[y_site-1,y_site+1,0,freshnc[n_i].variables[var_nc].shape[3]])[0][:,1,:]
    z_r_nutri = depths.get_zr_zw_tind(nutrinc[n_i],grid_nc,0,[y_site-1,y_site+1,0,nutrinc[n_i].variables[var_nc].shape[3]])[0][:,1,:]
    z_r_fulll = depths.get_zr_zw_tind(fulllnc[n_i],grid_nc,0,[y_site-1,y_site+1,0,fulllnc[n_i].variables[var_nc].shape[3]])[0][:,1,:]
    
    z_r_cntrl[z_r_cntrl>1E10] = np.nan
    z_r_fresh[z_r_fresh>1E10] = np.nan
    z_r_nutri[z_r_nutri>1E10] = np.nan
    z_r_fulll[z_r_fulll>1E10] = np.nan
    
    #ax.plot(-1*h_slice,color='k') # bottom depth  
    p_plot_cntrl = ax.flat[0].pcolor(lon_slice,z_r_cntrl,roms_var_cntrl,cmap=c_map,vmin=v_min,vmax=v_max)
    p_plot_fresh = ax.flat[1].pcolor(lon_slice,z_r_fresh,roms_var_fresh,cmap=c_map,vmin=v_min,vmax=v_max)
    p_plot_nutri = ax.flat[2].pcolor(lon_slice,z_r_nutri,roms_var_nutri,cmap=c_map,vmin=v_min,vmax=v_max)
    p_plot_fulll = ax.flat[3].pcolor(lon_slice,z_r_fulll,roms_var_fulll,cmap=c_map,vmin=v_min,vmax=v_max)
    
    p0 = ax.flat[1].get_position().get_points().flatten()
    p1 = ax.flat[3].get_position().get_points().flatten()
    cb_ax = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])
    
    cb = fig.colorbar(p_plot_cntrl,cax=cb_ax,orientation='vertical',format='%.1f')
    cb.set_label(cblabel,fontsize=axis_tick_size)
    cb.ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)
    # mark pipe location
    if loc == 'OCSD' or loc == 'JWPCP':
        ax.flat[0].scatter(lon_site+.006,z_r_cntrl[0,p_loc],marker='^',facecolors='orange',s=psize)
        ax.flat[1].scatter(lon_site+.006,z_r_fresh[0,p_loc],marker='^',facecolors='orange',s=psize)
        ax.flat[2].scatter(lon_site+.006,z_r_nutri[0,p_loc],marker='^',facecolors='orange',s=psize)
        ax.flat[3].scatter(lon_site+.006,z_r_fulll[0,p_loc],marker='^',facecolors='orange',s=psize)
    if loc == 'HTP' or loc == 'PLWTP':
        ax.flat[0].scatter(lon_site,dpipe,marker='^',facecolors='orange',s=psize)
        ax.flat[1].scatter(lon_site,dpipe,marker='^',facecolors='orange',s=psize)
        ax.flat[2].scatter(lon_site,dpipe,marker='^',facecolors='orange',s=psize)
        ax.flat[3].scatter(lon_site,dpipe,marker='^',facecolors='orange',s=psize)
    ax.flat[0].set_xlim([ind_st,ind_en])
    ax.flat[1].set_xlim([ind_st,ind_en])
    ax.flat[2].set_xlim([ind_st,ind_en])
    ax.flat[3].set_xlim([ind_st,ind_en])
    ax.flat[0].set_ylim([dp_st,dp_en])
    ax.flat[1].set_ylim([dp_st,dp_en])
    ax.flat[2].set_ylim([dp_st,dp_en])
    ax.flat[3].set_ylim([dp_st,dp_en])
    
    # plot contours
    clinecolor = 'k'
    c_plt_cntrl = ax.flat[0].contour(lon_reshape[:,ind_st_p:ind_en_p],z_r_cntrl[:,ind_st_p:ind_en_p],roms_var_cntrl[:,ind_st_p:ind_en_p],clines,colors=clinecolor,linewidths=1)
    c_plt_fresh = ax.flat[1].contour(lon_reshape[:,ind_st_p:ind_en_p],z_r_fresh[:,ind_st_p:ind_en_p],roms_var_fresh[:,ind_st_p:ind_en_p],clines,colors=clinecolor,linewidths=1)
    c_plt_nutri = ax.flat[2].contour(lon_reshape[:,ind_st_p:ind_en_p],z_r_nutri[:,ind_st_p:ind_en_p],roms_var_nutri[:,ind_st_p:ind_en_p],clines,colors=clinecolor,linewidths=1)
    c_plt_fulll = ax.flat[3].contour(lon_reshape[:,ind_st_p:ind_en_p],z_r_fulll[:,ind_st_p:ind_en_p],roms_var_fulll[:,ind_st_p:ind_en_p],clines,colors=clinecolor,linewidths=1)
    
    ax.flat[0].clabel(c_plt_cntrl,fontsize=9,fmt='%.1f',inline=1)
    ax.flat[1].clabel(c_plt_fresh,fontsize=9,fmt='%.1f',inline=1)
    ax.flat[2].clabel(c_plt_nutri,fontsize=9,fmt='%.1f',inline=1)
    ax.flat[3].clabel(c_plt_fulll,fontsize=9,fmt='%.1f',inline=1)
    
    
    ax.flat[0].set_ylabel('Depth (m)',fontsize=axis_tick_size)
    ax.flat[2].set_ylabel('Depth (m)',fontsize=axis_tick_size)
    
    ax.flat[2].set_xlabel('Longitude',fontsize=axis_tick_size)
    ax.flat[3].set_xlabel('Longitude',fontsize=axis_tick_size)
    
    # tick spacing 
    tick_spacingx = 0.04
    ax.flat[0].xaxis.set_major_locator(ticker.MultipleLocator(tick_spacingx))
    ax.flat[1].xaxis.set_major_locator(ticker.MultipleLocator(tick_spacingx))
    ax.flat[2].xaxis.set_major_locator(ticker.MultipleLocator(tick_spacingx))
    ax.flat[3].xaxis.set_major_locator(ticker.MultipleLocator(tick_spacingx))
    
    tick_spacingy = 15
    ax.flat[0].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacingy))
    ax.flat[1].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacingy))
    ax.flat[2].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacingy))
    ax.flat[3].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacingy))
    
    # remove labels
    ax.flat[1].get_yaxis().set_ticklabels([])
    ax.flat[3].get_yaxis().set_ticklabels([])
    
    ax.flat[0].get_xaxis().set_ticklabels([])
    ax.flat[1].get_xaxis().set_ticklabels([])
    
    
    fig.suptitle(calendar.month_name[int(savename[n_i][savename[n_i].index('M')+1:savename[n_i].index('M')+1+2])]+' '+savename[n_i][savename[n_i].index('Y')+1:savename[n_i].index('Y')+1+4]+' Average '+loc,fontsize=axis_tick_size)
    
    ax.flat[0].set_title('CTRL',fontsize=axis_tick_size)
    ax.flat[1].set_title('FRESH',fontsize=axis_tick_size)
    ax.flat[2].set_title('NUTR',fontsize=axis_tick_size)
    ax.flat[3].set_title('FULL',fontsize=axis_tick_size)
    
    ax.flat[0].tick_params(axis='both',which='major',labelsize=axis_tick_size)
    ax.flat[1].tick_params(axis='both',which='major',labelsize=axis_tick_size)
    ax.flat[2].tick_params(axis='both',which='major',labelsize=axis_tick_size)
    ax.flat[3].tick_params(axis='both',which='major',labelsize=axis_tick_size)
    
    
    fig.savefig(savepath+savename[n_i],bbox_inches='tight')
    print(savename)
    plt.close()

