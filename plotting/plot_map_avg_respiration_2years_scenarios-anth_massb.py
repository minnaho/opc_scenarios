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

#plt.ion()

timep = 'fullts'
# contour % line
clinemin = -5
clinemax = 5
clinemin2 = -10
clinemax2 = 10

savepath = './figs/2years/'

exp = ['cntrl_initap_realistic',
       'fulll_2012_2017',
       'PNDN_only_realistic',
       'pndn50_fixriver',
       'pndn90_fixriver',
       'FNDN_only_realistic',
       'fndn50_fixriver',
       'fndn90_fixriver']

title_exp = ['CTRL',
             'ANTH',
             '50% N Red.',
             '50% N Red.\n50% Recy.',
             '50% N Red.\n90% Recy.',
             '85% N Red.',
             '85% N Red.\n50% Recy.',
             '85% N Red.\n90% Recy.']

# roms var
varstr = 'respiration'
cblabel = '$\Delta$ '+varstr+' mmol m$^{-2}$'

s2d = 86400

# color map
#c_map = cmocean.cm.dense
#c_map = cmocean.cm.thermal
#c_map = cmocean.cm.algae
#c_map = cmocean.cm.deep
#c_map = cmocean.cm.delta
#c_map = 'PRGn'
c_map = cmocean.cm.balance
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
v_max = 10
v_min = -10

###########################
# mass balance
###########################
massbpath = '/data/project3/minnaho/opc_scenarios/plotting/massbalance_L2/newmb/budget_ww/'
varn = 'O2'
matn = 'MATBGCF'
matc = 'MATVARC'
matp = 'MATPHYF'

fst = 'outputs_'+varn+'_'
dep = '_0-200m'

# anth
matrn = h5py.File(massbpath+fst+'L2SCB_AP'+dep+'/'+matn+'.mat','r') 
matrp = h5py.File(massbpath+fst+'L2SCB_AP'+dep+'/'+matp+'.mat','r') 
datan = matrn.get(matn)
datap = matrp.get(matp)

# get dates
datemat = np.squeeze(h5py.File(massbpath+fst+'L2SCB_AP'+dep+'/'+matc+'.mat','r').get(matc)['date'])
dt = pd.to_datetime(datemat-719529, unit='D')

# get bgc variables O2
loss = np.squeeze(datan['LOSS'])
graze = np.squeeze(datan['GRAZE'])
remin = np.squeeze(datan['REMIN'])
sedre = np.squeeze(datan['SED_REMIN'])
ammox = np.squeeze(datan['AMMOX'])
nit = np.squeeze(datan['NIT'])

# get phys variables
adx = np.squeeze(datap['ADX'])
ady = np.squeeze(datap['ADY'])
adz = np.squeeze(datap['ADZ_bot'])+np.squeeze(datap['ADZ_top'])

# calculate respiration
respir_calc = loss+graze+remin+sedre+ammox+nit
print('anth 2016',str(np.nanmean(respir_calc[:12])*s2d))
print('anth 2017',str(np.nanmean(respir_calc[12:])*s2d))

# calculate phys terms
phys_tot = adx+ady+adz

# fullts
if timep == 'fullts':
    bgc_tot_anth = np.nanmean(respir_calc,axis=0)
    phys_tot_anth = np.nanmean(phys_tot,axis=0)
# spring
if timep == 'spring':
    bgc_tot_anth = np.nanmean((np.concatenate((respir_calc[5:8],respir_calc[17:20]))),axis=0)
    phys_tot_anth = np.nanmean((np.concatenate((phys_tot[5:8],phys_tot[17:20]))),axis=0)
# summer
if timep == 'summer':
    bgc_tot_anth = np.nanmean((np.concatenate((respir_calc[8:11],respir_calc[20:23]))),axis=0)
    phys_tot_anth = np.nanmean((np.concatenate((phys_tot[8:11],phys_tot[20:23]))),axis=0)
# winter
if timep == 'winter':
    bgc_tot_anth = np.nanmean((np.concatenate((respir_calc[2:5],respir_calc[14:17]))),axis=0)
    phys_tot_anth = np.nanmean((np.concatenate((phys_tot[2:5],phys_tot[14:17]))),axis=0)

# cntrl
matrn = h5py.File(massbpath+fst+'cntrl_initap_realistic'+dep+'/'+matn+'.mat','r') 
matrp = h5py.File(massbpath+fst+'cntrl_initap_realistic'+dep+'/'+matp+'.mat','r') 
datan = matrn.get(matn)
datap = matrp.get(matp)

# get dates
datemat = np.squeeze(h5py.File(massbpath+fst+'L2SCB_AP'+dep+'/'+matc+'.mat','r').get(matc)['date'])
dt = pd.to_datetime(datemat-719529, unit='D')

# get bgc variables O2
loss = np.squeeze(datan['LOSS'])
graze = np.squeeze(datan['GRAZE'])
remin = np.squeeze(datan['REMIN'])
sedre = np.squeeze(datan['SED_REMIN'])
ammox = np.squeeze(datan['AMMOX'])
nit = np.squeeze(datan['NIT'])

# get phys variables
adx = np.squeeze(datap['ADX'])
ady = np.squeeze(datap['ADY'])
adz = np.squeeze(datap['ADZ_bot'])+np.squeeze(datap['ADZ_top'])

# calculate respiration
respir_calc = loss+graze+remin+sedre+ammox+nit
print('cntrl 2016',str(np.nanmean(respir_calc[:12])*s2d))
print('cntrl 2017',str(np.nanmean(respir_calc[12:])*s2d))

# calculate phys terms
phys_tot = adx+ady+adz

# fullts
if timep == 'fullts':
    bgc_tot_cntrl = np.nanmean(respir_calc,axis=0)
    phys_tot_cntrl = np.nanmean(phys_tot,axis=0)
# spring
if timep == 'spring':
    bgc_tot_cntrl = np.nanmean((np.concatenate((respir_calc[5:8],respir_calc[17:20]))),axis=0)
    phys_tot_cntrl = np.nanmean((np.concatenate((phys_tot[5:8],phys_tot[17:20]))),axis=0)
# summer
if timep == 'summer':
    bgc_tot_cntrl = np.nanmean((np.concatenate((respir_calc[8:11],respir_calc[20:23]))),axis=0)
    phys_tot_cntrl = np.nanmean((np.concatenate((phys_tot[8:11],phys_tot[20:23]))),axis=0)
# winter
if timep == 'winter':
    bgc_tot_cntrl = np.nanmean((np.concatenate((respir_calc[2:5],respir_calc[14:17]))),axis=0)
    phys_tot_cntrl = np.nanmean((np.concatenate((phys_tot[2:5],phys_tot[14:17]))),axis=0)

npp_pos_mask = np.ones((len(exp)-2,mask_nc.shape[0],mask_nc.shape[1]))*np.nan
npp_neg_mask = np.ones((len(exp)-2,mask_nc.shape[0],mask_nc.shape[1]))*np.nan

# average of onshore vs offshore value
onshore_var = np.ones((len(exp)-2))*np.nan
offshore_var = np.ones((len(exp)-2))*np.nan

# 25th/75th percentile or 10/90 or 5/95
onshore_std_lower = np.ones((len(exp)-2))*np.nan
offshore_std_lower = np.ones((len(exp)-2))*np.nan
onshore_std_upper = np.ones((len(exp)-2))*np.nan
offshore_std_upper = np.ones((len(exp)-2))*np.nan

# average of onshore vs offshore positive and negative values 
onshore_var_pos = np.ones((len(exp)-2))*np.nan
offshore_var_pos = np.ones((len(exp)-2))*np.nan

onshore_var_neg = np.ones((len(exp)-2))*np.nan
offshore_var_neg = np.ones((len(exp)-2))*np.nan

# alongshore average
alongshore_var = np.ones((len(exp)-2,mask_nc.shape[1]))

# area in km^2 of positive and negative values
pos_area_var = np.ones((len(exp)-2))*np.nan
neg_area_var = np.ones((len(exp)-2))*np.nan

# area in km^2 onshore/offshore and positive/negative
pos_area_onshore = np.ones((len(exp)-2))*np.nan
neg_area_onshore = np.ones((len(exp)-2))*np.nan
pos_area_offshore =np.ones((len(exp)-2))*np.nan
neg_area_offshore =np.ones((len(exp)-2))*np.nan

# plot anth-cntrl, NM-cntrl
fig,ax = plt.subplots(1,3,figsize=[figw-2,4.1],sharey=True,subplot_kw=dict(projection=ccrs.PlateCarree()))
ax_i = 0
for e_i in [1,2,5]:
    # read data
    if e_i == 1:
        matrn = h5py.File(massbpath+fst+'L2SCB_AP'+dep+'/'+matn+'.mat','r')
        matrp = h5py.File(massbpath+fst+'L2SCB_AP'+dep+'/'+matp+'.mat','r')
    else:
        matrn = h5py.File(massbpath+fst+exp[e_i]+dep+'/'+matn+'.mat','r')
        matrp = h5py.File(massbpath+fst+exp[e_i]+dep+'/'+matp+'.mat','r')
        datemat = np.squeeze(h5py.File(massbpath+fst+exp[e_i]+dep+'/'+matc+'.mat','r').get(matc)['date'])

    datan = matrn.get(matn)
    datap = matrp.get(matp)

    # get dates

    # get bgc variables N
    loss = np.squeeze(datan['LOSS'])
    graze = np.squeeze(datan['GRAZE'])
    remin = np.squeeze(datan['REMIN'])
    sedre = np.squeeze(datan['SED_REMIN'])
    ammox = np.squeeze(datan['AMMOX'])
    nit = np.squeeze(datan['NIT'])
    # calculate bgc terms
    respir_calc = loss+graze+remin+sedre+ammox+nit

    bgc_tot_calc = respir_calc

    if timep == 'fullts':
        bgc_tot_sce = np.nanmean(respir_calc,axis=0)
        phys_tot_sce = np.nanmean(phys_tot,axis=0)

    bgc_tot_diff = (bgc_tot_sce - bgc_tot_cntrl)*s2d

    print('maps',title_exp[e_i])
    varpltbgc = bgc_tot_diff

    p_plot1 = ax.flat[ax_i].pcolormesh(lon_nc,lat_nc,varpltbgc,transform=ccrs.PlateCarree(),cmap=c_map,norm=mcolors.TwoSlopeNorm(vmin=-15,vcenter=0,vmax=15))
    #p_plot1 = ax.flat[ax_i].pcolormesh(lon_nc,lat_nc,varpltbgc,transform=ccrs.PlateCarree(),cmap=c_map,norm=mcolors.TwoSlopeNorm(vcenter=0))
    #c_plot = ax.flat[ax_i].contour(lon_nc,lat_nc,varpltbgc,[clinemax],transform=ccrs.PlateCarree(),colors='red',linestyles='dashed')
    c_plot = ax.flat[ax_i].contour(lon_nc,lat_nc,varpltbgc,[clinemax2],transform=ccrs.PlateCarree(),colors='red',linestyles='solid')
    #c_plot = ax.flat[ax_i].contour(lon_nc,lat_nc,varpltbgc,[clinemin],transform=ccrs.PlateCarree(),colors='blue',linestyles='dashed')
    c_plot = ax.flat[ax_i].contour(lon_nc,lat_nc,varpltbgc,[clinemin2],transform=ccrs.PlateCarree(),colors='blue',linestyles='solid')

    ax.flat[ax_i].set_title(title_exp[e_i],fontsize=axfont)
    ax_i += 1

for a_i in range(len(ax.flat)):
    # mark pipe location
    #for l_i in range(len(lon_potw)):
    #    ax.flat[a_i].scatter(lon_potw[l_i],lat_potw[l_i],marker='o',facecolors='none',edgecolors='blue',s=100)
    # other grid stuff
    ax.flat[a_i].add_feature(coast_10m,facecolor='None',edgecolor='k')
    ax.flat[a_i].add_feature(cpf.BORDERS,facecolor='None',edgecolor='k')
    ax.flat[a_i].set_extent(extent)
    gl = ax.flat[a_i].gridlines(draw_labels={'bottom':'x','left':'y'},xlabel_style=dict(size=axfont),ylabel_style=dict(size=axfont))
    if a_i >= 1:
        gl.left_labels=False


# colorbar
p0 = ax.flat[2].get_position().get_points().flatten()
p1 = ax.flat[2].get_position().get_points().flatten()
cb_ax = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])

cb1 = fig.colorbar(p_plot1,cax=cb_ax,orientation='vertical',ticks=mticker.MultipleLocator(5))
cb1.set_label(cblabel,fontsize=axfont)
cb1.ax.tick_params(axis='both',which='major',labelsize=axfont)

savename1 = 'map_2years-anth_100m_bgc_NM_'+varstr+'_'+region_name+'_abs_'+timep+'.png'

fig.savefig(savepath+savename1,bbox_inches='tight')

# plot each scenario-anth

fig1,ax1 = plt.subplots(2,3,figsize=[figw,figh],subplot_kw=dict(projection=ccrs.PlateCarree()))
fig2,ax2 = plt.subplots(2,3,figsize=[figw,figh],subplot_kw=dict(projection=ccrs.PlateCarree()))

for e_i in range(2,len(exp)):
    # read data
    matrn = h5py.File(massbpath+fst+exp[e_i]+dep+'/'+matn+'.mat','r') 
    matrp = h5py.File(massbpath+fst+exp[e_i]+dep+'/'+matp+'.mat','r') 
    datan = matrn.get(matn)
    datap = matrp.get(matp)
    # get dates
    datemat = np.squeeze(h5py.File(massbpath+fst+exp[e_i]+dep+'/'+matc+'.mat','r').get(matc)['date'])
    # get bgc variables N
    loss = np.squeeze(datan['LOSS'])
    graze = np.squeeze(datan['GRAZE'])
    remin = np.squeeze(datan['REMIN'])
    sedre = np.squeeze(datan['SED_REMIN'])
    ammox = np.squeeze(datan['AMMOX'])
    nit = np.squeeze(datan['NIT'])
    # get phys variables
    adx = np.squeeze(datap['ADX'])
    ady = np.squeeze(datap['ADY'])
    adz = np.squeeze(datap['ADZ_bot'])+np.squeeze(datap['ADZ_top'])
    # calculate bgc terms
    respir_calc = loss+graze+remin+sedre+ammox+nit
    print(exp[e_i],'2016',str(np.nanmean(respir_calc[:12])*s2d))
    print(exp[e_i],'2017',str(np.nanmean(respir_calc[12:])*s2d))
 
    bgc_tot_calc = respir_calc
    
    # calculate phys terms
    phys_tot = adx+ady+adz

    # fullts
    if timep == 'fullts':
        bgc_tot_sce = np.nanmean(respir_calc,axis=0)
        phys_tot_sce = np.nanmean(phys_tot,axis=0)
    # spring
    if timep == 'spring':
        bgc_tot_sce = np.nanmean((np.concatenate((respir_calc[5:8],respir_calc[17:20]))),axis=0)
        phys_tot_sce = np.nanmean((np.concatenate((phys_tot[5:8],phys_tot[17:20]))),axis=0)
    # summer
    if timep == 'summer':
        bgc_tot_sce = np.nanmean((np.concatenate((respir_calc[8:11],respir_calc[20:23]))),axis=0)
        phys_tot_sce = np.nanmean((np.concatenate((phys_tot[8:11],phys_tot[20:23]))),axis=0)
    # winter
    if timep == 'winter':
        bgc_tot_sce = np.nanmean((np.concatenate((respir_calc[2:5],respir_calc[14:17]))),axis=0)
        phys_tot_sce = np.nanmean((np.concatenate((phys_tot[2:5],phys_tot[14:17]))),axis=0)

    bgc_tot_diff = (bgc_tot_sce - bgc_tot_anth)*s2d
    phys_tot_diff = (phys_tot_sce - phys_tot_anth)*s2d
    
    bgc_tot_pos = bgc_tot_diff*npp_pos_mask[e_i-2]
    bgc_tot_neg = bgc_tot_diff*npp_neg_mask[e_i-2]
    
    bgc_sum_pos = np.nansum(bgc_tot_pos)
    bgc_sum_neg = np.nansum(bgc_tot_neg)
    
    phys_tot_pos = phys_tot_diff*npp_pos_mask[e_i-2]
    phys_tot_neg = phys_tot_diff*npp_neg_mask[e_i-2]
    
    phys_sum_pos = np.nansum(phys_tot_pos)
    phys_sum_neg = np.nansum(phys_tot_neg)
    
    print('maps',title_exp[e_i])
    varpltbgc = bgc_tot_diff
    varpltphys = phys_tot_diff
    #onshore_var[e_i-2] = data_onshore - anth_onshore
    #offshore_var[e_i-2] = data_offshore - anth_offshore

    # interannual variability
    # difference of scenario from cntrl
    #expcdiff = datanc-cntrld
    # difference of scenario from anth-cntrl
    #expddiff = ((expcdiff - diffmonth)/diffmonth)*100
    

    # calculate area of change > 5 and < -5 mmol m-2 d-1
    # in km2
    #npp_pos = np.where(varplt>=clinemax)
    #npp_neg = np.where(varplt<=clinemin)
    #pos_area = np.nansum(sizex[npp_pos[0],npp_pos[1]]*sizey[npp_pos[0],npp_pos[1]])
    #neg_area = np.nansum(sizex[npp_neg[0],npp_neg[1]]*sizey[npp_neg[0],npp_neg[1]])
    #print('area >'+str(clinemax)+' change NPP: '+str(int(pos_area))+' km^2')
    #print('area <'+str(clinemin)+' change NPP: '+str(int(neg_area))+' km^2')
    #pos_area_var[e_i-2] = pos_area
    #neg_area_var[e_i-2] = neg_area
    #npp_pos_mask[e_i-2,npp_pos[0],npp_pos[1]] = 1
    #npp_neg_mask[e_i-2,npp_neg[0],npp_neg[1]] = 1

    ## get anth values of where scenario > anth and scenario < anth
    #anth_onshore_pos = np.nanmean((anthavg*mask_onshore)[npp_pos]) 
    #anth_onshore_neg = np.nanmean((anthavg*mask_onshore)[npp_neg]) 
    #anth_offshore_pos = np.nanmean((anthavg*mask_offshore)[npp_pos]) 
    #anth_offshore_neg = np.nanmean((anthavg*mask_offshore)[npp_neg]) 

    ## calculate mean of onshore vs offshore where scenario > anth
    ## and scenario < anth
    #onshore_var_pos[e_i-2] = np.nanmean((dataavg*mask_onshore)[npp_pos]) - anth_onshore_pos
    #onshore_var_neg[e_i-2] = np.nanmean((dataavg*mask_onshore)[npp_neg]) - anth_onshore_neg
    #offshore_var_pos[e_i-2] = np.nanmean((dataavg*mask_offshore)[npp_pos]) - anth_offshore_pos
    #offshore_var_neg[e_i-2] = np.nanmean((dataavg*mask_offshore)[npp_neg]) - anth_offshore_neg
    #

    ## alongshore profile
    #alongshore_var[e_i-2] = np.nanmean(((dataavg*mask_onshore)-(anthavg*mask_onshore)),axis=0)

    p_plot1 = ax1.flat[e_i-2].pcolormesh(lon_nc,lat_nc,varpltbgc,transform=ccrs.PlateCarree(),cmap=c_map,norm=mcolors.TwoSlopeNorm(vmin=v_min,vcenter=0,vmax=v_max))
    #p_plot1 = ax1.flat[e_i-2].pcolormesh(lon_nc,lat_nc,varpltbgc,transform=ccrs.PlateCarree(),cmap=c_map,norm=mcolors.TwoSlopeNorm(vcenter=0))
    c_plot = ax1.flat[e_i-2].contour(lon_nc,lat_nc,varpltbgc,[clinemax],transform=ccrs.PlateCarree(),colors='red',linestyles='dashed')
    c_plot = ax1.flat[e_i-2].contour(lon_nc,lat_nc,varpltbgc,[clinemax2],transform=ccrs.PlateCarree(),colors='red',linestyles='solid')
    c_plot = ax1.flat[e_i-2].contour(lon_nc,lat_nc,varpltbgc,[clinemin],transform=ccrs.PlateCarree(),colors='blue',linestyles='dashed')
    c_plot = ax1.flat[e_i-2].contour(lon_nc,lat_nc,varpltbgc,[clinemin2],transform=ccrs.PlateCarree(),colors='blue',linestyles='solid')

    p_plot2 = ax2.flat[e_i-2].pcolormesh(lon_nc,lat_nc,varpltphys,transform=ccrs.PlateCarree(),cmap=c_map,norm=mcolors.TwoSlopeNorm(vmin=v_min,vcenter=0,vmax=v_max))
    #p_plot2 = ax2.flat[e_i-2].pcolormesh(lon_nc,lat_nc,varpltphys,transform=ccrs.PlateCarree(),cmap=c_map,norm=mcolors.TwoSlopeNorm(vcenter=0))
    c_plot = ax2.flat[e_i-2].contour(lon_nc,lat_nc,varpltphys,[clinemax],transform=ccrs.PlateCarree(),colors='lightgreen')
    c_plot = ax2.flat[e_i-2].contour(lon_nc,lat_nc,varpltphys,[clinemax2],transform=ccrs.PlateCarree(),colors='darkgreen')
    c_plot = ax2.flat[e_i-2].contour(lon_nc,lat_nc,varpltphys,[clinemin],transform=ccrs.PlateCarree(),colors='mediumpurple')
    c_plot = ax2.flat[e_i-2].contour(lon_nc,lat_nc,varpltphys,[clinemin2],transform=ccrs.PlateCarree(),colors='purple')
    
    ax1.flat[e_i-2].set_title(title_exp[e_i],fontsize=axfont)
    ax2.flat[e_i-2].set_title(title_exp[e_i],fontsize=axfont)

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

    ax2.flat[a_i].add_feature(coast_10m,facecolor='None',edgecolor='k')
    ax2.flat[a_i].add_feature(cpf.BORDERS,facecolor='None',edgecolor='k')
    ax2.flat[a_i].set_extent(extent)
    gl = ax2.flat[a_i].gridlines(draw_labels={'bottom':'x','left':'y'},xlabel_style=dict(size=axfont),ylabel_style=dict(size=axfont))
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

p0 = ax2.flat[2].get_position().get_points().flatten()
p1 = ax2.flat[5].get_position().get_points().flatten()
cb_ax2 = fig2.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])

cb2 = fig2.colorbar(p_plot2,cax=cb_ax2,orientation='vertical')
cb2.set_label(cblabel,fontsize=axfont)
cb2.ax.tick_params(axis='both',which='major',labelsize=axfont)

#savename = 'map_2years_100m_int_'+varstr+'_'+region_name+'_'+exp[-1]+'_abs_fullts.png'
savename1 = 'map_2years-anth_100m_bgc_'+varstr+'_'+region_name+'_abs_'+timep+'.png'
savename2 = 'map_2years-anth_100m_phys_'+varstr+'_'+region_name+'_abs_'+timep+'.png'

fig1.savefig(savepath+savename1,bbox_inches='tight')
fig2.savefig(savepath+savename2,bbox_inches='tight')

'''
###########################
# onshore vs offshore bars
###########################
print('bars')
fig,ax = plt.subplots(2,1,figsize=[12,5])
ax.flat[0].bar(range(len(title_exp[2:])),onshore_var)
ax.flat[1].bar(title_exp[2:],offshore_var)
fig.supylabel(varstr+' mmol m$^{-2}$',fontsize=axfont)
ax.flat[0].set_xticks([])
fig.suptitle('Scenario - ANTH',fontsize=axfont)
ax.flat[0].set_title('Onshore',fontsize=axfont)
ax.flat[1].set_title('Offshore',fontsize=axfont)
ax.flat[0].tick_params(axis='both',which='major',labelsize=axfont)
ax.flat[1].tick_params(axis='both',which='major',labelsize=axfont)
fig.tight_layout()
fig.savefig(savepath+'bar_'+varstr+'_onshore_offshore_'+timep+'.png',bbox_inches='tight')

fig,ax = plt.subplots(4,1,figsize=[12,8])
ax.flat[0].bar(range(len(title_exp[2:])),onshore_var_pos)
ax.flat[1].bar(range(len(title_exp[2:])),onshore_var_neg)
ax.flat[2].bar(title_exp[2:],offshore_var_pos)
ax.flat[3].bar(title_exp[2:],offshore_var_neg)
ax.flat[0].set_xticks([])
ax.flat[1].set_xticks([])
ax.flat[2].set_xticks([])
fig.suptitle('Scenario - ANTH',fontsize=axfont)
ax.flat[0].set_title('Onshore Positive',fontsize=axfont)
ax.flat[1].set_title('Onshore Negative',fontsize=axfont)
ax.flat[2].set_title('Offshore Positive',fontsize=axfont)
ax.flat[3].set_title('Offshore Negative',fontsize=axfont)
ax.flat[0].tick_params(axis='both',which='major',labelsize=axfont)
ax.flat[1].tick_params(axis='both',which='major',labelsize=axfont)
ax.flat[2].tick_params(axis='both',which='major',labelsize=axfont)
ax.flat[3].tick_params(axis='both',which='major',labelsize=axfont)
fig.supylabel(varstr+' mmol m$^{-2}$',fontsize=axfont)
fig.tight_layout()
fig.savefig(savepath+'bar_pos_neg_'+varstr+'_onshore_offshore_'+timep+'.png',bbox_inches='tight')

############################################
# total km npp > clinemax and npp < clinemin
############################################
fig,ax = plt.subplots(2,1,figsize=[12,8])
ax.flat[0].bar(range(len(title_exp[2:])),pos_area_var)
ax.flat[1].bar(title_exp[2:],neg_area_var)
fig.supylabel('km',fontsize=axfont)
ax.flat[0].set_xticks([])
fig.suptitle('Scenario - ANTH',fontsize=axfont)
ax.flat[0].set_title('Area of Increased NPP > 5 mmol m$^{-2}$',fontsize=axfont)
ax.flat[1].set_title('Area of Decreased NPP < -5 mmol m$^{-2}$',fontsize=axfont)
ax.flat[0].tick_params(axis='both',which='major',labelsize=axfont)
ax.flat[1].tick_params(axis='both',which='major',labelsize=axfont)
fig.tight_layout()
fig.savefig(savepath+'bar_'+varstr+'_area_'+timep+'.png',bbox_inches='tight')

####################
# alongshore profile
###################
latavg = np.nanmean(lat_nc*mask_onshore,axis=0)
fig,ax = plt.subplots(1,1,figsize=[6,12])
for a_i in range(len(alongshore_var)):
    ax.plot(alongshore_var[a_i][:-30],range(len(alongshore_var[a_i][:-30])))
'''
