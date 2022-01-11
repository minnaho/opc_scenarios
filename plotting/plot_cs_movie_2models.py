################################################
# plot cross section of particles
# with residence time x days
# consider all y points (not just one latitude line))
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
var_nc = 'w'
#var_nc = 'biomass'
#cblabel = 'Density (kg m$^{-3}$)'
cblabel = var_nc+' (m s$^{-1}$)'
#cblabel = 'Temperature C'
#cblabel = 'mmol C m$^{-3}$'
#cblabel = var_nc+' (mmol m$^{-3}$)'
#cblabel = var_nc+' (PSU)'

# color map
#c_map = cmocean.cm.dense
c_map = cmocean.cm.balance
#c_map = cmocean.cm.algae

# choose year and month
start_year = 1997
end_year = 2000

# between 1 and 12
start_month = 2
end_month = 12

# path of outputs
cntrlpath = '/data/project6/ROMS/L2SCB_1997_2000/DAILY/'
fulllpath = '/data/project6/ROMS/L2SCB_AP/daily/'

if var_nc == 'rho':
    rho0 = 1027.4
    clines = [1023.5,1024,1024.5,1025,1025.5,1026,1026.5]
if var_nc == 'NH4' or var_nc == 'NO3':
    clines = [0.5,1,5,10,13]
if var_nc == 'salt':
    clines = [32.9,33,33.1,33.2,33.3,33.4,33.5,33.6,33.7,33.8,33.9,34,34.1]
if var_nc == 'temp':
    clines = [9,10,11,12,13,14,15,16,17,18,19,20]
if var_nc == 'biomass':
    clines = [1,3,5,7,10,15,20]
if var_nc == 'w':
    clines = np.array([-0.1,-0.08,-0.06,-0.04,-0.02,.02,.04,.06,.08,.1])*1E-2


# roms file name
ncfile = []
savename = []

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
            ncfile.append('l2_scb_avg.'+year_month+'.nc')
            savename.append('cs_CTRL_ANTH_'+loc+'_'+var_nc+'_'+year_month+'.png')

# outputs
cntrlnc = []
fulllnc = []

for n_i in range(len(ncfile)):
    cntrlnc.append(cntrlpath+ncfile[n_i])
    fulllnc.append(fulllpath+ncfile[n_i])


# grid path
grid_nc = l2grid.grid_nc
lat_nc = l2grid.lat_nc
lon_nc = l2grid.lon_nc
h_nc = l2grid.h_nc

if loc == 'HTP':
    # north pipe
    lat_site = 33.920625 # pipe coordinates as implemented in model
    lon_site = -118.529712
    ind_st = -118.7 # range of lon to plot over
    ind_en = -118.47
    
    dp_st = -200 # set depth of domain
    dp_en = 0

    dpipe = -60 # depth of pipe
    
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
if loc == 'OCSan':
    lat_site = 33.57667
    lon_site = -118.01
    ind_st = -118.1
    ind_en = -117.963
    
    dp_st = -70
    dp_en = 0
    
    dpipe = -55 

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
lon_slice_l = list(lon_slice)*Dataset(cntrlnc[0],'r').variables['temp'].shape[1]
lon_reshape = np.array(lon_slice_l).reshape(Dataset(cntrlnc[0],'r').variables['temp'].shape[1],lon_nc.shape[1]) 

# bounds to draw contours
ind_st_p = np.nanmin(np.unique(np.where((lon_reshape[:,:]>ind_st)&(lon_reshape<ind_en))[1]))
ind_en_p = np.nanmax(np.unique(np.where((lon_reshape[:,:]>ind_st)&(lon_reshape<ind_en))[1]))

figw = 14
figh = 7.5

msize = 50
psize = 300

axis_tick_size = 14

for n_i in range(len(ncfile)):
    fig,ax = plt.subplots(1,2,figsize=[figw,figh])
    # roms field from output at y_site slice
    if var_nc == 'rho':
        roms_var_cntrl = np.array(Dataset(cntrlnc[n_i],'r').variables[var_nc][0,:,y_site,:])+rho0
        roms_var_fulll = np.array(Dataset(fulllnc[n_i],'r').variables[var_nc][0,:,y_site,:])+rho0

    elif var_nc == 'biomass':
        roms_var_cntrl_dtc = np.array(Dataset(cntrlnc[n_i],'r').variables['DIATC'][0,:,y_site,:])
        roms_var_cntrl_dtc[roms_var_cntrl_dtc>1E10] = 0
        roms_var_cntrl_spc = np.array(Dataset(cntrlnc[n_i],'r').variables['SPC'][0,:,y_site,:])
        roms_var_cntrl_spc[roms_var_cntrl_spc>1E10] = 0
        roms_var_cntrl_dzc = np.array(Dataset(cntrlnc[n_i],'r').variables['DIAZC'][0,:,y_site,:])
        roms_var_cntrl_dzc[roms_var_cntrl_dzc>1E10] = 0

        roms_var_nutri_dzc[roms_var_nutri_dzc>1E10] = 0

        roms_var_fulll_dtc = np.array(Dataset(fulllnc[n_i],'r').variables['DIATC'][0,:,y_site,:])
        roms_var_fulll_dtc[roms_var_fulll_dtc>1E10] = 0
        roms_var_fulll_spc = np.array(Dataset(fulllnc[n_i],'r').variables['SPC'][0,:,y_site,:])
        roms_var_fulll_spc[roms_var_fulll_spc>1E10] = 0
        roms_var_fulll_dzc = np.array(Dataset(fulllnc[n_i],'r').variables['DIAZC'][0,:,y_site,:])
        roms_var_fulll_dzc[roms_var_fulll_dzc>1E10] = 0

        roms_var_cntrl = roms_var_cntrl_dtc+roms_var_cntrl_spc+roms_var_cntrl_dzc
        roms_var_fresh = roms_var_fresh_dtc+roms_var_fresh_spc+roms_var_fresh_dzc
        roms_var_nutri = roms_var_nutri_dtc+roms_var_nutri_spc+roms_var_nutri_dzc
        roms_var_fulll = roms_var_fulll_dtc+roms_var_fulll_spc+roms_var_fulll_dzc

    else:
        roms_var_cntrl = np.array(Dataset(cntrlnc[n_i],'r').variables[var_nc][0,:,y_site,:])
        roms_var_fulll = np.array(Dataset(fulllnc[n_i],'r').variables[var_nc][0,:,y_site,:])
    
    roms_var_cntrl[roms_var_cntrl>1E10] = np.nan
    roms_var_fulll[roms_var_fulll>1E10] = np.nan
    
    # max and min of color bar
    if var_nc == 'rho':
        v_max = 1026.5
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
    if var_nc == 'biomass':
        v_max = 15
        #v_min = np.nanmin(roms_var_cntrl)
        v_min = 0
    if var_nc == 'w':
        v_max = 1E-3
        v_min = -1E-3

    
    
    # get depths at this y_site slice for each scenario
    z_r_cntrl = depths.get_zr_zw_tind(Dataset(cntrlnc[n_i],'r'),grid_nc,0,[y_site-1,y_site+1,0,Dataset(cntrlnc[n_i],'r').variables['temp'].shape[3]])[0][:,1,:]
    z_r_fulll = depths.get_zr_zw_tind(Dataset(fulllnc[n_i],'r'),grid_nc,0,[y_site-1,y_site+1,0,Dataset(fulllnc[n_i],'r').variables['temp'].shape[3]])[0][:,1,:]
    
    z_r_cntrl[z_r_cntrl>1E10] = np.nan
    z_r_fulll[z_r_fulll>1E10] = np.nan
    
    #ax.plot(-1*h_slice,color='k') # bottom depth  
    p_plot_cntrl = ax.flat[0].pcolor(lon_slice,z_r_cntrl,roms_var_cntrl,cmap=c_map,vmin=v_min,vmax=v_max)
    p_plot_fulll = ax.flat[1].pcolor(lon_slice,z_r_fulll,roms_var_fulll,cmap=c_map,vmin=v_min,vmax=v_max)
    
    p0 = ax.flat[1].get_position().get_points().flatten()
    p1 = ax.flat[1].get_position().get_points().flatten()
    cb_ax = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])
    
    if var_nc != 'w':
        cb = fig.colorbar(p_plot_cntrl,cax=cb_ax,orientation='vertical',format='%.1f')
    else:
        cb = fig.colorbar(p_plot_cntrl,cax=cb_ax,orientation='vertical',format='%.0e')
    cb.set_label(cblabel,fontsize=axis_tick_size)
    cb.ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)
    # mark pipe location
    if loc == 'OCSan' or loc == 'JWPCP':
        ax.flat[0].scatter(lon_site+.006,dpipe,marker='^',facecolors='orange',s=psize)
        ax.flat[1].scatter(lon_site+.006,dpipe,marker='^',facecolors='orange',s=psize)
    if loc == 'HTP' or loc == 'PLWTP':
        ax.flat[0].scatter(lon_site,dpipe,marker='^',facecolors='orange',s=psize)
        ax.flat[1].scatter(lon_site,dpipe,marker='^',facecolors='orange',s=psize)
    ax.flat[0].set_xlim([ind_st,ind_en])
    ax.flat[1].set_xlim([ind_st,ind_en])
    ax.flat[0].set_ylim([dp_st,dp_en])
    ax.flat[1].set_ylim([dp_st,dp_en])
    
    # plot contours
    fmtc = ticker.LogFormatterSciNotation()
    fmtc.create_dummy_axis()
    clinecolor = 'k'
    c_plt_cntrl = ax.flat[0].contour(lon_reshape[:,ind_st_p:ind_en_p],z_r_cntrl[:,ind_st_p:ind_en_p],roms_var_cntrl[:,ind_st_p:ind_en_p],clines,colors=clinecolor,linewidths=1,fmt=fmtc)
    c_plt_fulll = ax.flat[1].contour(lon_reshape[:,ind_st_p:ind_en_p],z_r_fulll[:,ind_st_p:ind_en_p],roms_var_fulll[:,ind_st_p:ind_en_p],clines,colors=clinecolor,linewidths=1,fmt=fmtc)
    
    if var_nc != 'w':
        ax.flat[0].clabel(c_plt_cntrl,fontsize=9,fmt='%.1f',inline=1)
        ax.flat[1].clabel(c_plt_fulll,fontsize=9,fmt='%.1f',inline=1)
    else:
        ax.flat[0].clabel(c_plt_cntrl,fontsize=9,fmt='%.0e',inline=1)
        ax.flat[1].clabel(c_plt_fulll,fontsize=9,fmt='%.0e',inline=1)
        
    
    ax.flat[0].set_ylabel('Depth (m)',fontsize=axis_tick_size)
    
    ax.flat[0].set_xlabel('Longitude',fontsize=axis_tick_size)
    ax.flat[1].set_xlabel('Longitude',fontsize=axis_tick_size)
    
    # tick spacing 
    tick_spacingx = 0.1
    ax.flat[0].xaxis.set_major_locator(ticker.MultipleLocator(tick_spacingx))
    ax.flat[1].xaxis.set_major_locator(ticker.MultipleLocator(tick_spacingx))
    
    tick_spacingy = 25
    ax.flat[0].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacingy))
    ax.flat[1].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacingy))
    
    # remove labels
    ax.flat[1].get_yaxis().set_ticklabels([])
    
    fig.suptitle(savename[n_i][savename[n_i].index('D')+1:savename[n_i].index('D')+1+2]+' '+calendar.month_name[int(savename[n_i][savename[n_i].index('M')+1:savename[n_i].index('M')+1+2])]+' '+savename[n_i][savename[n_i].index('Y')+1:savename[n_i].index('Y')+1+4]+' Average '+loc,fontsize=axis_tick_size)
    
    ax.flat[0].set_title('CTRL',fontsize=axis_tick_size)
    ax.flat[1].set_title('ANTH',fontsize=axis_tick_size)
    
    ax.flat[0].tick_params(axis='both',which='major',labelsize=axis_tick_size)
    ax.flat[1].tick_params(axis='both',which='major',labelsize=axis_tick_size)
    
    fig.savefig(savepath+savename[n_i],bbox_inches='tight')
    print(savename[n_i])
    plt.close()

