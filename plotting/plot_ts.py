# plot timeseries of integrated and
# sliced varibles
import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import l2grid
import numpy as np
from netCDF4 import Dataset,num2date
import matplotlib.pyplot as plt

plt.ion()

var = 'O2'
dp = '50m'
tp = 'sli'
perc = True

if dp == '50m':
    ncpath = '../postprocessing/ts_int_sli/'
if dp == 'srf':
    ncpath = '../postprocessing/surf_sli/'

figpath = './figs/ts/'

if perc == True:
    ylabel = '% '+var+' change'
    # biomass
    #y_lim0 = -75
    #y_lim1 = 180
    # O2
    y_lim0 = -20
    y_lim1 = 20
if perc == False:
    ylabel = 'mmol change'

savename = 'ts_change_'+var+'_'+dp+'.png'


fresh_units = 'days since 1999-07-01'

# control NH4, NO3,O2,biomass start from 1997-02-01
# control carbonate (pH,omegas) start from 1999-07-01

# masks
mask_path = '/data/project1/minnaho/make_masks/mask_scb.nc'
mask_smm = np.array(Dataset(mask_path,'r').variables['mask_sm'])

fresh = ncpath+tp+'_'+dp+'_fresh_'+var+'.nc'
nutri = ncpath+tp+'_'+dp+'_nutri_'+var+'.nc'
cntrl = ncpath+tp+'_'+dp+'_cntrl_'+var+'.nc'

#if var == 'O2' or var == 'omega_co2sys' or var == 'omega_juranek' or var == 'pH':
#    depth_type = 'slice'
#    fresh = ncpath+'sli_'+dp+'_fresh_'+var+'.nc'
#    nutri = ncpath+'sli_'+dp+'_nutri_'+var+'.nc'
#    cntrl = ncpath+'sli_'+dp+'_cntrl_'+var+'.nc'
#else:
#    depth_type = 'integrated'
#    fresh = ncpath+'int_'+dp+'_fresh_'+var+'.nc'
#    nutri = ncpath+'int_'+dp+'_nutri_'+var+'.nc'
#    cntrl = ncpath+'int_'+dp+'_cntrl_'+var+'.nc'

if var == 'omega_co2sys' or var == 'omega_juranek' or var == 'pH':
    freshnc = np.array(Dataset(fresh,'r').variables[var])
    nutrinc = np.array(Dataset(nutri,'r').variables[var])
    cntrlnc_long = np.array(Dataset(cntrl,'r').variables[var])
    cntrl_units = fresh_units
if var == 'totbiomass':
    freshnc = np.squeeze(np.array(Dataset(fresh,'r').variables['totc']))
    nutrinc = np.squeeze(np.array(Dataset(nutri,'r').variables['totc']))
    cntrlnc_long = np.squeeze(np.array(Dataset(cntrl,'r').variables['totc']))
    cntrl_units = fresh_units

else:
    freshnc = np.array(Dataset(fresh,'r').variables['var'])
    nutrinc = np.array(Dataset(nutri,'r').variables['var'])
    cntrlnc_long = np.array(Dataset(cntrl,'r').variables['var'])
    cntrl_units = 'days since 1997-02-01'

freshnc[freshnc>1E10] = np.nan
nutrinc[nutrinc>1E10] = np.nan
cntrlnc_long[cntrlnc_long>1E10] = np.nan

# date - subsample cntrl because not same length
# Y1999M07 - Y2000M10
fresh_dt = num2date(np.arange(freshnc.shape[0]),fresh_units,only_use_cftime_datetimes=False,only_use_python_datetimes=True)

cntrl_dt = num2date(np.arange(cntrlnc_long.shape[0]),cntrl_units,only_use_cftime_datetimes=False,only_use_python_datetimes=True)

# find time periods 
cntrlnc = cntrlnc_long[np.where((cntrl_dt>=fresh_dt[0])&(cntrl_dt<=fresh_dt[-1]))]

# multiply by mask then average over region
freshmask = np.nanmean(np.nanmean(freshnc*mask_smm,axis=1),axis=1)
nutrimask = np.nanmean(np.nanmean(nutrinc*mask_smm,axis=1),axis=1)
cntrlmask = np.nanmean(np.nanmean(cntrlnc*mask_smm,axis=1),axis=1)

# calculate % change
if perc == True:
    freshch = ((freshmask-cntrlmask)/cntrlmask)*100
    nutrich = ((nutrimask-cntrlmask)/cntrlmask)*100
if perc == False:
    freshch = freshmask-cntrlmask
    nutrich = nutrimask-cntrlmask

# plotting
axisfont = 16
figw = 14
figh = 10

fig,ax = plt.subplots(2,1,figsize=[figw,figh])
ax[0].plot(fresh_dt,freshch,color='k')
ax[0].fill_between(fresh_dt,freshch,where=freshch<0,color='blue',alpha=0.5)
ax[0].fill_between(fresh_dt,freshch,where=freshch>0,color='yellow',alpha=0.5)
ax[0].set_title(var+' '+tp+' '+dp+' FRESH-CNTRL',fontsize=axisfont)
ax[0].plot(fresh_dt,np.ones(fresh_dt.shape[0])*0,color='k',linestyle='--')

ax[1].plot(fresh_dt,nutrich,color='k')
ax[1].plot(fresh_dt,np.ones(fresh_dt.shape[0])*0,color='k',linestyle='--')
ax[1].fill_between(fresh_dt,nutrich,where=nutrich<0,color='blue',alpha=0.5)
ax[1].fill_between(fresh_dt,nutrich,where=nutrich>0,color='yellow',alpha=0.5)
ax[1].set_title(var+' '+tp+' '+dp+' NUTRI-CNTRL',fontsize=axisfont)

for i in range(len(ax)):
    ax[i].set_ylabel(ylabel,fontsize=axisfont)
    ax[i].tick_params(axis='both',which='major',labelsize=axis_tick_font)
    ax[i].yaxis.set_ticks_position('both')
    ax[i].xaxis.set_ticks_position('both')
    if perc == True:
        ax[i].set_ylim([y_lim0,y_lim1])


fig.savefig(figpath+savename,bbox_inches='tight')
