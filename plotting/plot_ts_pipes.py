# plot timeseries of integrated and
# sliced varibles
# includes pipes only run
import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import l2grid
import numpy as np
from netCDF4 import Dataset,num2date
import matplotlib.pyplot as plt

plt.ion()

var = 'biomass'
dp = '100m'
tp = 'int'
msk = 'mask_sp'
#msk = 'mask_sm'
perc = True

if dp == '50m' or dp == '100m':
    ncpath = '../postprocessing/ts_int_sli/'
if dp == 'srf':
    ncpath = '../postprocessing/surf_sli/'

figpath = './figs/ts/'

if perc == True:
    savename = 'ts_change_'+var+'_'+dp+'_'+msk[-2:]+'_perc.png'
    ylabel = '% '+var+' change'
    # biomass
    y_lim0 = -75
    y_lim1 = 180
    # O2
    #y_lim0 = -25
    #y_lim1 = 25
if perc == False:
    savename = 'ts_change_'+var+'_'+dp+'_'+msk[-2:]+'_abso.png'
    ylabel = 'mmol change'



fresh_units = 'days since 1999-07-01'

# control NH4, NO3,O2,biomass start from 1997-02-01
# control carbonate (pH,omegas) start from 1999-07-01

# masks
mask_path = '/data/project1/minnaho/make_masks/mask_scb.nc'
mask_smm = np.array(Dataset(mask_path,'r').variables[msk])

fresh = ncpath+tp+'_'+dp+'_fresh_'+var+'.nc'
nutri = ncpath+tp+'_'+dp+'_nutri_'+var+'.nc'
cntrl = ncpath+tp+'_'+dp+'_cntrl_'+var+'.nc'
pipes = ncpath+tp+'_'+dp+'_pipes_'+var+'.nc'
fulll = ncpath+tp+'_'+dp+'_fulll_'+var+'.nc'

if var == 'omega_co2sys' or var == 'omega_juranek' or var == 'pH':
    freshnc = np.array(Dataset(fresh,'r').variables[var])
    nutrinc = np.array(Dataset(nutri,'r').variables[var])
    pipesnc = np.array(Dataset(pipes,'r').variables[var])
    fulllnc = np.array(Dataset(fulll,'r').variables[var])
    cntrlnc_long = np.array(Dataset(cntrl,'r').variables[var])
    cntrl_units = fresh_units
if var == 'totbiomass':
    freshnc = np.squeeze(np.array(Dataset(fresh,'r').variables['totc']))
    nutrinc = np.squeeze(np.array(Dataset(nutri,'r').variables['totc']))
    pipesnc = np.squeeze(np.array(Dataset(pipes,'r').variables['totc']))
    fulllnc = np.squeeze(np.array(Dataset(fulll,'r').variables['totc']))
    cntrlnc_long = np.squeeze(np.array(Dataset(cntrl,'r').variables['totc']))
    cntrl_units = fresh_units
if var == 'biomass':
    freshnc = np.array(Dataset(fresh,'r').variables['var'])
    nutrinc = np.array(Dataset(nutri,'r').variables['var'])
    pipesnc = np.array(Dataset(pipes,'r').variables['var'])
    fulllnc = np.array(Dataset(fulll,'r').variables['var'])
    cntrlnc_long = np.array(Dataset(cntrl,'r').variables['var'])
    cntrl_units = fresh_units
else:
    freshnc = np.array(Dataset(fresh,'r').variables['var'])
    nutrinc = np.array(Dataset(nutri,'r').variables['var'])
    pipesnc = np.array(Dataset(pipes,'r').variables['var'])
    fulllnc = np.array(Dataset(fulll,'r').variables['var'])
    cntrlnc_long = np.array(Dataset(cntrl,'r').variables['var'])
    cntrl_units = 'days since 1997-02-01'

freshnc[freshnc>1E10] = np.nan
nutrinc[nutrinc>1E10] = np.nan
pipesnc[pipesnc>1E10] = np.nan
fulllnc[fulllnc>1E10] = np.nan
cntrlnc_long[cntrlnc_long>1E10] = np.nan

# date - subsample cntrl because not same length
# Y1999M07 - Y2000M09
# remove last 31 days from Oct to plot pipes
plt_dt = num2date(np.arange(pipesnc.shape[0]),fresh_units,only_use_cftime_datetimes=False,only_use_python_datetimes=True)

cntrl_dt = num2date(np.arange(cntrlnc_long.shape[0]),cntrl_units,only_use_cftime_datetimes=False,only_use_python_datetimes=True)

# find time periods 
freshnc = freshnc[:plt_dt.shape[0]]
nutrinc = nutrinc[:plt_dt.shape[0]]
pipesnc = pipesnc[:plt_dt.shape[0]]
fulllnc = fulllnc[:plt_dt.shape[0]]
cntrlnc = cntrlnc_long[np.where((cntrl_dt>=plt_dt[0])&(cntrl_dt<=plt_dt[-1]))]

# multiply by mask then average over region
freshmask = np.nanmean(np.nanmean(freshnc*mask_smm,axis=1),axis=1)
nutrimask = np.nanmean(np.nanmean(nutrinc*mask_smm,axis=1),axis=1)
pipesmask = np.nanmean(np.nanmean(pipesnc*mask_smm,axis=1),axis=1)
fulllmask = np.nanmean(np.nanmean(fulllnc*mask_smm,axis=1),axis=1)
cntrlmask = np.nanmean(np.nanmean(cntrlnc*mask_smm,axis=1),axis=1)

# calculate % change
if perc == True:
    freshch = ((freshmask-cntrlmask)/cntrlmask)*100
    nutrich = ((nutrimask-cntrlmask)/cntrlmask)*100
    pipesch = ((pipesmask-cntrlmask)/cntrlmask)*100
    fulllch = ((fulllmask-cntrlmask)/cntrlmask)*100
if perc == False:
    freshch = freshmask-cntrlmask
    nutrich = nutrimask-cntrlmask
    pipesch = pipesmask-cntrlmask
    fulllch = fulllmask-cntrlmask

# plotting
axisfont = 16
figw = 14
figh = 16

fig,ax = plt.subplots(4,1,figsize=[figw,figh])
ax[0].plot(plt_dt,freshch,color='k')
ax[0].fill_between(plt_dt,freshch,where=freshch<0,color='blue',alpha=0.5)
ax[0].fill_between(plt_dt,freshch,where=freshch>0,color='yellow',alpha=0.5)
ax[0].set_title(var+' '+tp+' '+dp+' FRESH-CNTRL',fontsize=axisfont)
ax[0].plot(plt_dt,np.ones(plt_dt.shape[0])*0,color='k',linestyle='--')

ax[1].plot(plt_dt,nutrich,color='k')
ax[1].plot(plt_dt,np.ones(plt_dt.shape[0])*0,color='k',linestyle='--')
ax[1].fill_between(plt_dt,nutrich,where=nutrich<0,color='blue',alpha=0.5)
ax[1].fill_between(plt_dt,nutrich,where=nutrich>0,color='yellow',alpha=0.5)
ax[1].set_title(var+' '+tp+' '+dp+' NUTRI-CNTRL',fontsize=axisfont)

ax[2].plot(plt_dt,pipesch,color='k')
ax[2].plot(plt_dt,np.ones(plt_dt.shape[0])*0,color='k',linestyle='--')
ax[2].fill_between(plt_dt,pipesch,where=pipesch<0,color='blue',alpha=0.5)
ax[2].fill_between(plt_dt,pipesch,where=pipesch>0,color='yellow',alpha=0.5)
ax[2].set_title(var+' '+tp+' '+dp+' PIPES-CNTRL',fontsize=axisfont)

ax[3].plot(plt_dt,fulllch,color='k')
ax[3].plot(plt_dt,np.ones(plt_dt.shape[0])*0,color='k',linestyle='--')
ax[3].fill_between(plt_dt,fulllch,where=fulllch<0,color='blue',alpha=0.5)
ax[3].fill_between(plt_dt,fulllch,where=fulllch>0,color='yellow',alpha=0.5)
ax[3].set_title(var+' '+tp+' '+dp+' FULL-CNTRL',fontsize=axisfont)

for i in range(len(ax)):
    ax[i].set_ylabel(ylabel,fontsize=axisfont)
    ax[i].tick_params(axis='both',which='major',labelsize=axisfont)
    ax[i].yaxis.set_ticks_position('both')
    ax[i].xaxis.set_ticks_position('both')
    ax[i].set_ylim([ax[0].get_ylim()[0],ax[1].get_ylim()[1]])
    #if perc == True:
    #    ax[i].set_ylim([y_lim0,y_lim1])

fig.tight_layout()

fig.savefig(figpath+savename,bbox_inches='tight')
