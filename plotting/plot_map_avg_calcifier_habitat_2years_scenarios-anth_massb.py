################################################
# plot maps of opc scenarios - 8 maps total
# 1 month or season
# 7 OPC scenarios
# 1 CTRL scenario
# conda activate cartopy_update 
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
import pandas as pd
import h5py
import scipy.io
import scipy.ndimage

#plt.ion()

timep = 'fullts'
# contour % line
clinemin = -10
clinemax = 10
clinemin2 = -20
clinemax2 = 20

savepath = './figs/2years/'

exp = ['cntrl_initap_realistic',
       'fulll_2012_2017',
       'PNDN_only_realistic',
       'pndn50_fixriver',
       'pndn90_fixriver',
       'FNDN_only_realistic',
       'fndn50_fixriver',
       'fndn90_fixriver']

exp1 = ['cntrl_initap',
       'loads1617',
       'PNDN_only',
       'pndn50',
       'pndn90',
       'FNDN_only',
       'fndn50',
       'fndn90']

title_exp = ['CTRL',
             'ANTH',
             '50% N Red.',
             '50% N Red.\n50% Recy.',
             '50% N Red.\n90% Recy.',
             '85% N Red.',
             '85% N Red.\n50% Recy.',
             '85% N Red.\n90% Recy.']

# roms var
varstr = 'Calcifier Habitat'
cblabel = '$\Delta$ '+varstr+' %'

s2d = 86400

# color map
#c_map = cmocean.cm.dense
#c_map = cmocean.cm.thermal
#c_map = cmocean.cm.algae
#c_map = cmocean.cm.deep
#c_map = cmocean.cm.delta
#c_map = 'PRGn'
c_map = cmocean.cm.balance_r
#c_map1 = cmocean.cm.algae
#c_map = 'PRGn'

# grid path
grid_nc = l2grid.grid_nc
lat_nc = l2grid.lat_nc
lon_nc = l2grid.lon_nc
h_nc = l2grid.h_nc
pm_nc = l2grid.pm_nc
pn_nc = l2grid.pn_nc
mask_nc = l2grid.mask_nc

# get size of cells
sizex = 1E-3/pm_nc
sizey = 1E-3/pn_nc


# mask name
region_mask = Dataset('/data/project1/minnaho/make_masks/mask_scb.nc','r')
mask_cst = np.array(region_mask.variables['mask_coast'])
mask_cst[mask_cst==0] = np.nan

mask_onshore = mask_cst
regtitle = '15 km Coast'

mask_temp = np.copy(mask_nc)
mask_temp[:,:20] = np.nan
mask_temp[:20,:] = np.nan
mask_temp[-20:,:] = np.nan
mask_grid = mask_temp
regtitle = 'Bightwide'

# do opposite of coastal band mask
regtitle = 'Offshore'
mask_temp = np.copy(mask_cst)
mask_temp[mask_temp==1] = 2
mask_temp[np.isnan(mask_temp)] = 1
mask_temp[mask_temp==2] = 0
mask_temp = mask_temp*mask_nc
mask_temp[:,:20] = np.nan
mask_temp[:20,:] = np.nan
mask_temp[-20:,:] = np.nan
mask_offshore = mask_temp
#regtitle = 'Offshore'
#mask7[mask9==1] = 1
#mask_mult = mask7

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

figw = 18
figh = 8

axfont = 16

# large pipes lat and lon
major_nc = Dataset('/data/project1/minnaho/potw_outfall_data/updated_2013_2017/major_potw_data/major_potw_1971_2017_monthly.nc','r')
lat_potw = np.array(major_nc.variables['latitude'])
lon_potw = np.array(major_nc.variables['longitude'])

coast_10m = cpf.NaturalEarthFeature('physical','coastline','10m')

# max and min of color bar
v_max = 20
v_min = -20

###########################
# read data
###########################
massbpath = 'calcifierhabitat_nutrientmngmt.mat'

matr = scipy.io.loadmat(massbpath)
# keys 
matkeys = ['Zt_pndn','Zt_p50', 'Zt_p90','Zt_fndn','Zt_f50', 'Zt_f90']   

# plot each scenario-anth

fig1,ax1 = plt.subplots(2,3,figsize=[figw,figh],subplot_kw=dict(projection=ccrs.PlateCarree()))

m_i = 0
for e_i in range(2,len(exp)):
    # read data
    respir_calc = matr[matkeys[m_i]] 
    varpltbgc = matr[matkeys[m_i]]

    p_plot1 = ax1.flat[e_i-2].pcolormesh(lon_nc,lat_nc,varpltbgc,transform=ccrs.PlateCarree(),cmap=c_map,norm=mcolors.TwoSlopeNorm(vmin=v_min,vcenter=0,vmax=v_max))
    #p_plot1 = ax1.flat[e_i-2].pcolormesh(lon_nc,lat_nc,varpltbgc,transform=ccrs.PlateCarree(),cmap=c_map,norm=mcolors.TwoSlopeNorm(vcenter=0))
    #varpltbgc = scipy.ndimage.zoom(varpltbgc,10)
    varpltbgc = scipy.ndimage.gaussian_filter(varpltbgc,sigma=2,order=0)
    c_plot = ax1.flat[e_i-2].contour(lon_nc,lat_nc,varpltbgc,[clinemax],transform=ccrs.PlateCarree(),colors='blue',linestyles='dashed')
    c_plot = ax1.flat[e_i-2].contour(lon_nc,lat_nc,varpltbgc,[clinemax2],transform=ccrs.PlateCarree(),colors='blue',linestyles='solid')
    c_plot = ax1.flat[e_i-2].contour(lon_nc,lat_nc,varpltbgc,[clinemin],transform=ccrs.PlateCarree(),colors='red',linestyles='dashed')
    c_plot = ax1.flat[e_i-2].contour(lon_nc,lat_nc,varpltbgc,[clinemin2],transform=ccrs.PlateCarree(),colors='red',linestyles='solid')
    c_h = ax1.flat[e_i-2].contour(lon_nc,lat_nc,h_nc,[200],transform=ccrs.PlateCarree(),colors='k',linestyles='solid')

    ax1.flat[e_i-2].set_title(title_exp[e_i],fontsize=axfont)
    m_i += 1

for a_i in range(len(ax1.flat)):
    # mark pipe location
    #for l_i in range(len(lon_potw)):
    #    ax.flat[a_i].scatter(lon_potw[l_i],lat_potw[l_i],marker='o',facecolors='none',edgecolors='blue',s=100)
    # other grid stuff
    ax1.flat[a_i].add_feature(coast_10m,facecolor='None',edgecolor='k')
    ax1.flat[a_i].add_feature(cpf.BORDERS,facecolor='None',edgecolor='k')
    ax1.flat[a_i].set_extent(extent)
    gl = ax1.flat[a_i].gridlines(draw_labels={'bottom':'x','left':'y'},xlabel_style=dict(size=axfont),ylabel_style=dict(size=axfont))
    if a_i >= 1 and a_i !=3 :
        gl.left_labels=False 
    if a_i >= 0 and a_i < 3 :
        gl.bottom_labels = False

# colorbar
p0 = ax1.flat[2].get_position().get_points().flatten()
p1 = ax1.flat[5].get_position().get_points().flatten()
cb_ax1 = fig1.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])

cb1 = fig1.colorbar(p_plot1,cax=cb_ax1,orientation='vertical')
cb1.set_label(cblabel,fontsize=axfont)
cb1.ax.tick_params(axis='both',which='major',labelsize=axfont)

savename1 = 'map_2years-anth_'+varstr+'_'+region_name+'_abs_'+timep+'.png'

fig1.savefig(savepath+savename1,bbox_inches='tight')

# plot histogram of area of calcifier habitat change
nbins = 100

cplt = ['orange','orange','orange','gray','gray','gray']
lsty = ['-','--',':','-','--',':']

title_exp_hist = ['50% N Red.',
             '50% N Red.\n50% Recy.',
             '50% N Red.\n90% Recy.',
             '85% N Red.',
             '85% N Red.\n50% Recy.',
             '85% N Red.\n90% Recy.']

fighist,axhist = plt.subplots(1,1,figsize=[8,6])

for e_i in range(3):
    varpltbgc = matr[matkeys[e_i]]
    n,bins = np.histogram(varpltbgc[~np.isnan(varpltbgc)],nbins)
    n = np.append(n,0)
    n = n*333*333/1E6
    axhist.plot(bins,n,color=cplt[e_i],linestyle=lsty[e_i],linewidth=3,label=title_exp_hist[e_i])

axhist.axvline(-10,color='red',linestyle='--')
axhist.axvline(10,color='blue',linestyle='--')
axhist.axvline(20,color='blue',linestyle='-')
axhist.set_ylim([0,4500])
axhist.set_xlim([-20,45])
axhist.set_ylabel('Area per % (km$^2$/%)',fontsize=axfont)
axhist.set_xlabel('Change in Calcifier Habitat (%)\n[relative to ANTH]',fontsize=axfont)
axhist.legend(loc='best',fontsize=axfont)
axhist.tick_params(axis='both',which='major',labelsize=axfont)
fighist.savefig(savepath+'hist_calcifier_habitat_pndn.png',bbox_inches='tight')

fighist,axhist = plt.subplots(1,1,figsize=[8,6])

for e_i in range(3,len(matkeys)):
    varpltbgc = matr[matkeys[e_i]]
    n,bins = np.histogram(varpltbgc[~np.isnan(varpltbgc)],nbins)
    n = np.append(n,0)
    n = n*333*333/1E6
    axhist.plot(bins,n,color=cplt[e_i],linestyle=lsty[e_i],linewidth=3,label=title_exp_hist[e_i])

axhist.axvline(-10,color='red',linestyle='--')
axhist.axvline(10,color='blue',linestyle='--')
axhist.axvline(20,color='blue',linestyle='-')
axhist.set_ylim([0,4500])
axhist.set_yticklabels([])
axhist.set_xlim([-20,45])
#axhist.set_ylabel('Area (km$^2$)',fontsize=axfont)
axhist.set_xlabel('Change in Calcifier Habitat (%)\n[relative to ANTH]',fontsize=axfont)
axhist.legend(loc='best',fontsize=axfont)
axhist.tick_params(axis='both',which='major',labelsize=axfont)
fighist.savefig(savepath+'hist_calcifier_habitat_fndn.png',bbox_inches='tight')
