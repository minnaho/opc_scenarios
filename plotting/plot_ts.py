# plot timeseries of integrated and
# sliced varibles
import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import l2grid
import numpy as np
from netCDF4 import Dataset,num2date
import matplotlib.pyplot as plt

#plt.ion()

region_name = 'coast'
var = 'O2'
var_nc = 'var'
dp = '50'
perc = True

ncpath = '/data/project6/minnaho/opc_scenarios/slices/'
if var == 'O2':
    filest = 'concat_sub_'+dp+'_'+var+'_'
else:
    filest = 'concat_'+dp+'_'+var+'_'

fileen = '.nc'

cntrlnc = 'concat_'+dp+'_'+var+'_cntrl'+fileen
cntrl_var = np.array(Dataset(ncpath+cntrlnc,'r').variables[var_nc])

figpath = './figs/ts/'

if perc == True:
    ylabel = '% '+var+' change'
    # biomass
    #y_lim0 = -75
    #y_lim1 = 180
    # O2
    y_lim0 = -150
    y_lim1 = 10
if perc == False:
    ylabel = 'mmol change'


# scenario names
exp = ['l1617','PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
title_exp = ['Loads 16-17','PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']

time_units = 'days since 1997-11-01'

dtplt = num2date(np.arange(cntrl_var.shape[0]),time_units,only_use_cftime_datetimes=False,only_use_python_datetimes=True)

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

mask_ssd[mask_ssd==0] = np.nan
mask_nsd[mask_nsd==0] = np.nan
mask_oc[mask_oc==0] = np.nan
mask_sp[mask_sp==0] = np.nan
mask_sm[mask_sm==0] = np.nan
mask_v[mask_v==0] = np.nan
mask_sb[mask_sb==0] = np.nan

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
if region_name == 'scb': # full L2 grid
    mask_mult = mask_nc
    regtitle = 'SCB'
if region_name == 'coast':
    mask_mult = mask_cst
    regtitle = '15 km coast'

# apply mask and average for control
cntrl_mask = cntrl_var*mask_mult
cntrl_avg = np.nanmean(np.nanmean(cntrl_mask,axis=1),axis=1)

# plotting
axisfont = 16
figw = 14
figh = 10

for n_i in range(len(exp)):

    roms_var = np.squeeze(np.array(Dataset(fpath[n_i],'r').variables[var_nc]))*mask_mult
    # multiply by mask then average over region
    roms_avg = np.nanmean(np.nanmean(roms_var,axis=1),axis=1)

    # calculate % change
    if perc == True:
        roms_plt = ((roms_avg-cntrl_avg)/cntrl_avg)*100
    if perc == False:
        roms_plt = roms_avg-cntrl_avg


    fig,ax = plt.subplots(1,1,figsize=[figw,figh])
    ax.plot(dtplt,roms_plt,color='k')
    ax.fill_between(dtplt,roms_plt,where=roms_plt<0,color='blue',alpha=0.5)
    ax.fill_between(dtplt,roms_plt,where=roms_plt>0,color='yellow',alpha=0.5)
    ax.set_title(var+' '+dp+'m '+title_exp[n_i],fontsize=axisfont)
    ax.plot(dtplt,np.ones(roms_plt.shape[0])*0,color='k',linestyle='--')

    ax.set_ylabel(ylabel,fontsize=axisfont)
    ax.tick_params(axis='both',which='major',labelsize=axisfont)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    if perc == True:
        ax.set_ylim([y_lim0,y_lim1])

    savename = 'ts_change_'+var+'_'+dp+'_'+region_name+'_'+exp[n_i]+'.png'
    fig.savefig(figpath+savename,bbox_inches='tight')
    plt.close()
