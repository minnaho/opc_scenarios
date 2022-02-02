# sum up all negative values in a mask
# used to see impact of treatment level on oxygen
# make my own 5-10 km mask circles?
import sys
import os
sys.path.append('/data/project3/minnaho/global/')
import numpy as np
from netCDF4 import Dataset,num2date
import matplotlib.pyplot as plt

# plot fresh vs nutrients vs control vs full
plt.ion()

savepath = './figs/scatter/'
region_name = 'coast'

# ROMS output location
outpath = '/data/project6/minnaho/opc_scenarios/ext_depth/'

# roms var
var_name = 'O2' 
var_nc = 'var' 
cblabel = 'mmol '+var_name+' m$^{-3}$'

year_month = 'summer1998'

# scenario names 
#exp = ['l1617','PNDN_only','FNDN_only','pndn50','pndn90','fndn50','fndn90']
#title_exp = ['Loads 16-17','PNDN only','FNDN only','PNDN 50','PNDN 90','FNDN 50','FNDN 90']
exp = ['l1617','PNDN_only','FNDN_only','pndn50','pndn90','fndn50','fndn90']
title_exp = ['Loads 16-17','PNDN only','FNDN only','PNDN 50','PNDN 90','FNDN 50','FNDN 90','OUTF']
#exp = ['l1617','PNDN_only','fndn90']
#title_exp = ['Loads 16-17','PNDN only','FNDN 90']

filest = 'sub_avg_'+year_month+'_0_80_'+var_name+'_'
fileen = '.nc'

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

fpath = []
for e_i in range(len(exp)):
    fpath.append(outpath+filest+exp[e_i]+fileen)

figw = 12
figh = 4

axis_tick_size = 14

savename = 'sum_neg_'+var_name+'_'+year_month+'_'+region_name+'_paige_nobounds.png'
fig,ax = plt.subplots(1,1,figsize=[figw,figh])

for n_i in range(len(exp)):

    roms_var = np.array(Dataset(fpath[n_i],'r').variables[var_nc])*mask_mult
    # sum up all negative values then make it a positive number to plot
    roms_neg = np.nansum(roms_var[roms_var<0])*-1

    ax.scatter(n_i,roms_neg)

# paige's run
paigerun = np.transpose(np.array(Dataset(outpath+'80m_outfall_negchange_O2_summer2000.nc','r').variables['O2_dif']),[0,2,1])
paigesum = np.nanmean(paigerun,axis=0)
paigemask = paigesum*mask_mult
# sum up all negative values then make it a positive number to plot
paigeplt = np.nansum(paigemask[paigemask<0])*-1

ax.scatter(n_i+1,paigeplt)

ax.set_yscale('log')
#ax.set_ybound(lower=3E5,upper=5E7)
ax.set_xticks(range(len(exp)+1))
ax.set_xticklabels(title_exp)
ax.set_xlabel('Scenario',fontsize=axis_tick_size)
ax.set_ylabel('Sum of O2 change',fontsize=axis_tick_size)
ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)

fig.suptitle('Sum of Negative O2 '+year_month+' Average '+regtitle,fontsize=axis_tick_size)


fig.savefig(savepath+savename,bbox_inches='tight')
print(savename)
#plt.close()

