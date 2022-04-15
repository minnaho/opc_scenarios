# SPATIAL AVERAGE
# sum up all negative values in a mask
# used to see impact of treatment level on oxygen
import sys
import os
sys.path.append('/data/project3/minnaho/global/')
import l2grid as l2grid
import numpy as np
from netCDF4 import Dataset,num2date
import matplotlib.pyplot as plt

# plot fresh vs nutrients vs control vs full
plt.ion()

savepath = './figs/scatter/'
region_name = 'grid'

# ROMS output location
outpath = '/data/project6/minnaho/opc_scenarios/ext_depth_200/'

# roms var
var_name = 'O2' 
var_nc = 'var' 
cblabel = 'mmol '+var_name+' m$^{-3}$'

year_month = 'spring1998'

filest = 'ext_0_200_'+var_name+'_'+year_month+'_'
fileen = '.nc'

# scenario names 
exp = ['l1617','PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
title_exp = ['Loads 16-17','PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']
#exp = ['l1617','PNDN_only','FNDN_only']
#title_exp = ['Loads 16-17','PNDN only','FNDN only']

# region masks
mask_nc = l2grid.mask_nc
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
if region_name == 'scb': # full L2 grid
    mask_mult = mask_nc
    regtitle = 'SCB'
if region_name == 'coast':
    mask_mult = mask_cst
    regtitle = '15 km coast'
if region_name == 'grid':
    mask_mult = mask_nc
    regtitle = 'full SCB'

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

fpath = []
for e_i in range(len(exp)):
    fpath.append(outpath+filest+exp[e_i]+fileen)

# read in control to subtract from
roms_cnt = np.squeeze(Dataset(outpath+filest+'cntrl'+fileen,'r').variables[var_nc])*mask_mult

figw = 12
figh = 4

axis_tick_size = 14

savename = 'sum_neg_'+var_name+'_'+year_month+'_'+region_name+'_200_'+exp[-1]+'_avg.png'
fig,ax = plt.subplots(1,1,figsize=[figw,figh])

for n_i in range(len(exp)):

    roms_var_read = np.squeeze(Dataset(fpath[n_i],'r').variables[var_nc])*mask_mult
    roms_var = roms_var_read - roms_cnt
    # sum up all negative values then make it a positive number to plot
    roms_var[roms_var>=0] = 0
    roms_avg = np.nanmean(np.nanmean(roms_var,axis=1),axis=1)
    roms_neg = np.nansum(roms_avg)*-1

    ax.bar(n_i,roms_neg)

#ax.set_yscale('log')
#ax.set_ybound(lower=3E5,upper=5E7)
ax.set_xticks(range(len(exp)))
ax.set_xticklabels(title_exp)
ax.set_xlabel('Scenario',fontsize=axis_tick_size)
ax.set_ylabel('Sum of O2 change mmol/m$^2$',fontsize=axis_tick_size)
ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)
ax.tick_params(axis='both',which='minor',labelsize=axis_tick_size)

fig.suptitle('Sum of Negative O2 '+year_month+' Average '+regtitle,fontsize=axis_tick_size)


fig.savefig(savepath+savename,bbox_inches='tight')
print(savename)
#plt.close()

