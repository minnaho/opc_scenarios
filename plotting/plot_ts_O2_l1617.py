# plot timeseries of integrated and
# sliced varibles
# compare to loads1617 scenario
import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import l2grid
import runmean
import numpy as np
from netCDF4 import Dataset,num2date
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

plt.ion()

region_name = 'mask3'
varstr = 'O2'
var_nc = 'var'
dp = '50'
perc = False

ncpath = '/data/project6/minnaho/opc_scenarios/slices/'
filest = 'concat_'+dp+'_'+varstr+'_'

fileen = '.nc'

cntrlnc = filest+'l1617'+fileen
cntrl_var = np.squeeze(Dataset(ncpath+cntrlnc,'r').variables[var_nc])

figpath = './figs/ts/'

if perc == True:
    ylabel = '% change'
else:
    ylabel = 'mmol change'
    

# scenario names
#exp = ['l1617','PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
#title_exp = ['Loads 16-17','PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']
exp = ['PNDN_only','pndn50','pndn90']
title_exp = ['PNDN only','PNDN 50','PNDN 90']
#exp = ['FNDN_only','fndn50','fndn90']
#title_exp = ['FNDN only','FNDN 50','FNDN 90']

time_units = 'days since 1997-11-01'

dtpltlong = num2date(np.arange(cntrl_var.shape[0]),time_units,only_use_cftime_datetimes=False,only_use_python_datetimes=True)

fpath = []
for e_i in range(len(exp)):
    fpath.append(ncpath+filest+exp[e_i]+fileen)

# region masks
region_mask = Dataset('/data/project1/minnaho/make_masks/mask_scb.nc','r')
mask_ssd = np.array(region_mask.variables['mask_ssd'])
mask_nsd = np.array(region_mask.variables['mask_nsd'])
mask_oc = np.array(region_mask.variables['mask_oc'])
mask_sp = np.array(region_mask.variables['mask_sp'])
mask_sm = np.array(region_mask.variables['mask_sm'])
mask_v = np.array(region_mask.variables['mask_v'])
mask_sb = np.array(region_mask.variables['mask_sb'])
# 15 km of coast
mask_cst = np.array(region_mask.variables['mask_coast'])

mask_nc = l2grid.mask_nc

mask_ssd[mask_ssd==0] = np.nan
mask_nsd[mask_nsd==0] = np.nan
mask_oc[mask_oc==0] = np.nan
mask_sp[mask_sp==0] = np.nan
mask_sm[mask_sm==0] = np.nan
mask_v[mask_v==0] = np.nan
mask_sb[mask_sb==0] = np.nan
mask_nc[mask_nc==0] = np.nan

# masks 0-9
masknum = Dataset('/data/project1/minnaho/make_masks/mask_gridL2.nc','r')
mask0 = np.array(masknum.variables['mask0'])
mask1 = np.array(masknum.variables['mask1'])
mask2 = np.array(masknum.variables['mask2'])
mask3 = np.array(masknum.variables['mask3'])
mask4 = np.array(masknum.variables['mask4'])
mask5 = np.array(masknum.variables['mask5'])
mask6 = np.array(masknum.variables['mask6'])
mask7 = np.array(masknum.variables['mask7'])
mask8 = np.array(masknum.variables['mask8'])
mask9 = np.array(masknum.variables['mask9'])

mask0[mask0==0] = np.nan
mask1[mask1==0] = np.nan
mask2[mask2==0] = np.nan
mask3[mask3==0] = np.nan
mask4[mask4==0] = np.nan
mask5[mask5==0] = np.nan
mask6[mask6==0] = np.nan
mask7[mask7==0] = np.nan
mask8[mask8==0] = np.nan
mask9[mask9==0] = np.nan


if region_name == 'ssd':
    mask_mult = mask_ssd
    regtitle = 'South San Diego'
if region_name == 'nsd':
    mask_mult = mask_nsd
    regtitle = 'North San Diego'
if region_name == 'oc':
    mask_mult = mask_oc
    regtitle = 'Orange County'
if region_name == 'sp':
    mask_mult = mask_sp
    regtitle = 'San Pedro'
if region_name == 'sm':
    mask_mult = mask_sm
    regtitle = 'Santa Monica Bay'
if region_name == 'v':
    mask_mult = mask_v
    regtitle = 'Ventura'
if region_name == 'sb':
    mask_mult = mask_sb
    regtitle = 'Santa Barbara'
if region_name == 'grid': # full L2 grid
    mask_mult = mask_nc
    regtitle = 'SCB'
if region_name == 'coast':
    mask_mult = mask_cst
    regtitle = '15 km coast'

if region_name == 'mask0':
    mask_mult = mask0
    regtitle = region_name
if region_name == 'mask1':
    mask_mult = mask1
    regtitle = region_name
if region_name == 'mask2':
    mask_mult = mask2
    regtitle = region_name
if region_name == 'mask3':
    mask_mult = mask3
    regtitle = region_name
if region_name == 'mask4':
    mask_mult = mask4
    regtitle = region_name
if region_name == 'mask5':
    mask_mult = mask5
    regtitle = region_name
if region_name == 'mask6':
    mask_mult = mask6
    regtitle = region_name
if region_name == 'mask7':
    mask_mult = mask7
    regtitle = region_name
if region_name == 'mask8':
    mask_mult = mask8
    regtitle = region_name
if region_name == 'mask9':
    mask_mult = mask9
    regtitle = region_name


# apply mask and average for control
cntrl_mask = cntrl_var*mask_mult
cntrl_avg = np.nanmean(np.nanmean(cntrl_mask,axis=1),axis=1)

# plotting
axisfont = 16
figw = 14
figh = 10

if perc == False:
    v_max = 25
    v_min = -25
else:
    v_max = 15
    v_min = -15

# moving average 7 days - shorten date array
dtplt = dtpltlong[1:-5]

fig,ax = plt.subplots(len(exp),1,figsize=[figw,figh])

for n_i in range(len(exp)):

    roms_var = np.squeeze(np.array(Dataset(fpath[n_i],'r').variables[var_nc]))*mask_mult
    # multiply by mask then average over region
    roms_avg = np.nanmean(np.nanmean(roms_var,axis=1),axis=1)

    if perc == True:
        roms_sub = ((roms_avg-cntrl_avg)/cntrl_avg)*100
    else:
        roms_sub = roms_avg-cntrl_avg

    roms_std = runmean.running_mean(np.nanstd(np.nanstd(roms_var,axis=1),axis=1),7)

    roms_plt = runmean.running_mean(roms_sub,7)

    ax[n_i].plot(dtplt,roms_plt,color='k')
    ax[n_i].plot(dtplt,roms_plt+roms_std,color='k',linestyle='--')
    ax[n_i].plot(dtplt,roms_plt-roms_std,color='k',linestyle='--')
    ax[n_i].fill_between(dtplt,roms_plt,where=roms_plt<0,color='blue',alpha=0.5)
    ax[n_i].fill_between(dtplt,roms_plt,where=roms_plt>0,color='yellow',alpha=0.5)
    ax[n_i].set_title(title_exp[n_i],fontsize=axisfont)
    ax[n_i].plot(dtplt,np.ones(roms_plt.shape[0])*0,color='k',linestyle='--')
    ax[n_i].set_ylim([v_min,v_max])
    #if perc == False:
    #    tick_spacingy = 100
    #else:
    #    tick_spacingy = 30
    #ax[n_i].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacingy))

    ax[n_i].set_ylabel(ylabel,fontsize=axisfont)
    ax[n_i].tick_params(axis='both',which='major',labelsize=axisfont)
    ax[n_i].yaxis.set_ticks_position('both')
    ax[n_i].xaxis.set_ticks_position('both')
    if n_i < len(exp)-1:
        ax[n_i].xaxis.set_ticklabels([])

if perc == True:
    savename = 'ts_change_'+varstr+'_'+dp+'_'+region_name+'_'+exp[n_i]+'_l1617_perc.png'
else:
    savename = 'ts_change_'+varstr+'_'+dp+'_'+region_name+'_'+exp[n_i]+'_l1617.png'

fig.savefig(figpath+savename,bbox_inches='tight')
