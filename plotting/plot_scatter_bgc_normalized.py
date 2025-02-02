# sum up all negative values in a mask
# used to see impact of treatment level on oxygen
# make my own 5-10 km mask circles?
import sys
import os
sys.path.append('/data/project3/minnaho/global/')
import l2grid
import numpy as np
from netCDF4 import Dataset,num2date
import glob as glob
import matplotlib.pyplot as plt

# plot fresh vs nutrients vs control vs full
plt.ion()

savepath = './figs/scatter/'
region_name = 'grid'

# ROMS output location
outpath = '/data/project6/minnaho/opc_scenarios/bgc_flux/'

# roms var
var_name = 'npp' 
var_nc = 'var_int' 
cblabel = 'mmol m$^{-2}$ d$^{-1}$'

#year_month = 'Y1998_M04_06'
year_month = 'fullts'
timename = '11/1997-11/1999'

# scenario names 
exp = ['PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
#title_exp = ['PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']
title_exp = ['50% N\nReduction','50% N\nReduction\n50% Recycle','50% N\nReduction\n90% Recycle','85% N\nReduction','85% N\nReduction\n50% Recycle','85% N\nReduction\n90% Recycle',]
#exp = ['PNDN_only','FNDN_only']
#title_exp = ['50% N\nReduction','85% N\nReduction']

#filest = 'int_avg_100m_50m_'
filest = 'avg_fullts_int_avg_100m_50m_'
#fileen = '_'+var_name+'_'+year_month+'.nc'
fileen = '_'+var_name+'.nc'

# region masks
mask_nc = l2grid.mask_nc
pm_nc = l2grid.pm_nc
pn_nc = l2grid.pn_nc

xisize = 1/pm_nc
etasize = 1/pn_nc

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

s2d = 86400

cntrl_var = np.squeeze(Dataset(outpath+filest+'cntrl'+fileen,'r').variables[var_nc])*mask_mult
cntrl_mean = np.nanmean(cntrl_var)*s2d
l1617_var = np.squeeze(Dataset(outpath+filest+'loads1617'+fileen,'r').variables[var_nc])*mask_mult
l1617_mean = np.nanmean(l1617_var)*s2d

potwr_var = np.squeeze(Dataset(outpath+filest+'POTW'+fileen,'r').variables[var_nc])*mask_mult
potwr_mean = np.nanmean(potwr_var)*s2d

# difference between l1617 and POTW only run = rivers 
potwd_var = np.squeeze(Dataset(outpath+filest+'loads1617-POTW'+fileen,'r').variables[var_nc])*mask_mult
potwd_mean = np.nanmean(potwd_var)*s2d

figw = 14
#figw = 10
figh = 4

axis_tick_size = 14

savename = var_name+'_'+year_month+'_'+region_name+'_'+exp[-1]+'_bar_norm.png'
fig,ax = plt.subplots(1,1,figsize=[figw,figh])

for n_i in range(len(exp)):

    roms_var = np.squeeze(Dataset(fpath[n_i],'r').variables[var_nc])*mask_mult
    roms_neg = np.nanmean(roms_var)*s2d

    print(exp[n_i],str(roms_neg))
    print(exp[n_i],'up to '+str(np.nanmin(((roms_var-l1617_var)/l1617_var)*100))+'% decrease in productivity')
    print(exp[n_i],'up to '+str(np.nanmax(((roms_var-l1617_var)/l1617_var)*100))+'% increase in productivity')
    ax.bar(n_i,roms_neg-cntrl_mean,bottom=cntrl_mean)

ax.set_ylim(bottom=cntrl_mean-(cntrl_mean*0.01),top=l1617_mean+(cntrl_mean*0.01))
ax.set_xticks(range(len(exp)))
ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*cntrl_mean,linestyle='--',color='purple')
ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*l1617_mean,linestyle='--',color='green')
#ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*potwd_mean,linestyle='--',color='blue')

ax.set_xticklabels(title_exp)
ax.set_xlabel('Scenario',fontsize=axis_tick_size)
ax.set_ylabel('Integrated NPP 100 m '+cblabel,fontsize=axis_tick_size)
ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)

#fig.suptitle('Integrated '+var_name+' '+timename+' '+regtitle,fontsize=axis_tick_size)


fig.savefig(savepath+savename,bbox_inches='tight')
print(savename)
#plt.close()

