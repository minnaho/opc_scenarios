import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import l2grid
import numpy as np
import scipy.stats as stats
from netCDF4 import Dataset,num2date,date2num
import matplotlib.pyplot as plt
from collections import Counter

roms_path = '/data/project6/minnaho/opc_scenarios/surf_sli/'
savepath = './figs/pdf/'

varnc = 'totchla'

varname = 'chla'

# nsd, ssd, oc, sp, sm, v, sb, or grid or coast (15 km)
region_name = 'coast'

savename = 'pdf_'+varname+'_'+region_name+'_log_cntrl_initap98'

# grid variables
mask_nc = l2grid.mask_nc


# read variables
fulll_nc1 = Dataset(roms_path+'concat_totchla_loads1617.nc','r')
fulll_nc2 = Dataset(roms_path+'concat_totchla_fulll_2012_2017.nc','r')
cntrl_nc1 = Dataset(roms_path+'concat_totchla_cntrl_initap.nc','r')
#cntrl_nc1 = Dataset(roms_path+'concat_totchla_cntrl.nc','r')
cntrl_nc2 = Dataset(roms_path+'concat_totchla_cntrl_2012_2017.nc','r')
pndnon_nc1 = Dataset(roms_path+'concat_totchla_PNDN_only.nc','r')
pndnon_nc2 = Dataset(roms_path+'concat_totchla_PNDN_only_realistic.nc','r')
fndnon_nc1 = Dataset(roms_path+'concat_totchla_FNDN_only.nc','r')
fndnon_nc2 = Dataset(roms_path+'concat_totchla_FNDN_only_realistic.nc','r')
pndn50_nc1 = Dataset(roms_path+'concat_totchla_pndn50.nc','r')
pndn50_nc2 = Dataset(roms_path+'concat_totchla_pndn50_realistic.nc','r')
pndn90_nc1 = Dataset(roms_path+'concat_totchla_pndn90.nc','r')
pndn90_nc2 = Dataset(roms_path+'concat_totchla_pndn90_realistic.nc','r')
fndn50_nc1 = Dataset(roms_path+'concat_totchla_fndn50.nc','r')
fndn50_nc2 = Dataset(roms_path+'concat_totchla_fndn50_realistic.nc','r')
fndn90_nc1 = Dataset(roms_path+'concat_totchla_fndn90.nc','r')
fndn90_nc2 = Dataset(roms_path+'concat_totchla_fndn90_realistic.nc','r')

# region masks
region_mask = Dataset('/data/project1/minnaho/make_masks/mask_scb.nc','r')
mask_ssd = np.array(region_mask.variables['mask_ssd'])
mask_nsd = np.array(region_mask.variables['mask_nsd'])
mask_oc = np.array(region_mask.variables['mask_oc'])
mask_sp = np.array(region_mask.variables['mask_sp'])
# my SM mask
mask_sm = np.array(region_mask.variables['mask_sm'])
# faycal's SM mask from PNAS paper
#mask_sm = np.transpose(np.array(h5py.File('../masksm.mat','r')['masksm']))
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
if region_name == 'grid': # full L2 grid
    mask_mult = mask_nc
    regtitle = 'SCB'
if region_name == 'coast':
    mask_mult = mask_cst
    regtitle = '15 km Coast'

fulll_mask1 = np.squeeze(fulll_nc1.variables[varnc])*mask_mult
fulll_mask2 = np.squeeze(fulll_nc2.variables[varnc])*mask_mult
cntrl_mask1 = np.squeeze(cntrl_nc1.variables[varnc])*mask_mult
cntrl_mask2 = np.squeeze(cntrl_nc2.variables[varnc])*mask_mult
pndnon_mask1 = np.squeeze(pndnon_nc1.variables[varnc])*mask_mult
pndnon_mask2 = np.squeeze(pndnon_nc2.variables[varnc])*mask_mult
fndnon_mask1 = np.squeeze(fndnon_nc1.variables[varnc])*mask_mult
fndnon_mask2 = np.squeeze(fndnon_nc2.variables[varnc])*mask_mult
pndn50_mask1 = np.squeeze(pndn50_nc1.variables[varnc])*mask_mult
pndn50_mask2 = np.squeeze(pndn50_nc2.variables[varnc])*mask_mult
pndn90_mask1 = np.squeeze(pndn90_nc1.variables[varnc])*mask_mult
pndn90_mask2 = np.squeeze(pndn90_nc2.variables[varnc])*mask_mult
fndn50_mask1 = np.squeeze(fndn50_nc1.variables[varnc])*mask_mult
fndn50_mask2 = np.squeeze(fndn50_nc2.variables[varnc])*mask_mult
fndn90_mask1 = np.squeeze(fndn90_nc1.variables[varnc])*mask_mult
fndn90_mask2 = np.squeeze(fndn90_nc2.variables[varnc])*mask_mult

fulll_mask1 = fulll_mask1[fulll_mask1>0]
fulll_mask2 = fulll_mask2[fulll_mask2>0]
cntrl_mask1 = cntrl_mask1[cntrl_mask1>0]
cntrl_mask2 = cntrl_mask2[cntrl_mask2>0]
pndnon_mask1 = pndnon_mask1[pndnon_mask1>0]
pndnon_mask2 = pndnon_mask2[pndnon_mask2>0]
fndnon_mask1 = fndnon_mask1[fndnon_mask1>0]
fndnon_mask2 = fndnon_mask2[fndnon_mask2>0]
pndn50_mask1 = pndn50_mask1[pndn50_mask1>0]
pndn50_mask2 = pndn50_mask2[pndn50_mask2>0]
pndn90_mask1 = pndn90_mask1[pndn90_mask1>0]
pndn90_mask2 = pndn90_mask2[pndn90_mask2>0]
fndn50_mask1 = fndn50_mask1[fndn50_mask1>0]
fndn50_mask2 = fndn50_mask2[fndn50_mask2>0]
fndn90_mask1 = fndn90_mask1[fndn90_mask1>0]
fndn90_mask2 = fndn90_mask2[fndn90_mask2>0]

fulll_var = np.concatenate((fulll_mask1,fulll_mask2))
cntrl_var = np.concatenate((cntrl_mask1,cntrl_mask2))
pndnon_var = np.concatenate((pndnon_mask1,pndnon_mask2))
fndnon_var = np.concatenate((fndnon_mask1,fndnon_mask2))
pndn50_var = np.concatenate((pndn50_mask1,pndn50_mask2))
pndn90_var = np.concatenate((pndn90_mask1,pndn90_mask2))
fndn50_var = np.concatenate((fndn50_mask1,fndn50_mask2))
fndn90_var = np.concatenate((fndn90_mask1,fndn90_mask2))

# remove all negative and 0 values

fulll_var[fulll_var>1E10] = np.nan
cntrl_var[cntrl_var>1E10] = np.nan
pndnon_var[pndnon_var>1E10] = np.nan
fndnon_var[fndnon_var>1E10] = np.nan
pndn50_var[pndn50_var>1E10] = np.nan
pndn90_var[pndn90_var>1E10] = np.nan
fndn50_var[fndn50_var>1E10] = np.nan
fndn90_var[fndn90_var>1E10] = np.nan

fulll_var[fulll_var<=0] = np.nan
cntrl_var[cntrl_var<=0] = np.nan
pndnon_var[pndnon_var<=0] = np.nan
fndnon_var[fndnon_var<=0] = np.nan
pndn50_var[pndn50_var<=0] = np.nan
pndn90_var[pndn90_var<=0] = np.nan
fndn50_var[fndn50_var<=0] = np.nan
fndn90_var[fndn90_var<=0] = np.nan

# stats
nbins = 500

n_fl,bins_fl,patch_fl = plt.hist(fulll_var.flatten(),nbins)
n_cl,bins_cl,patch_cl = plt.hist(cntrl_var.flatten(),nbins)
n_po,bins_po,patch_po = plt.hist(pndnon_var.flatten(),nbins)
n_fo,bins_fo,patch_fo = plt.hist(fndnon_var.flatten(),nbins)
n_p5,bins_p5,patch_p5 = plt.hist(pndn50_var.flatten(),nbins)
n_p9,bins_p9,patch_p9 = plt.hist(pndn90_var.flatten(),nbins)
n_f5,bins_f5,patch_f5 = plt.hist(fndn50_var.flatten(),nbins)
n_f9,bins_f9,patch_f9 = plt.hist(fndn90_var.flatten(),nbins)

# non nan values
flt_fl = np.where(~np.isnan(fulll_var.flatten()))[0].shape[0]
flt_cl = np.where(~np.isnan(cntrl_var.flatten()))[0].shape[0]
flt_po = np.where(~np.isnan(pndnon_var.flatten()))[0].shape[0]
flt_fo = np.where(~np.isnan(fndnon_var.flatten()))[0].shape[0]
flt_p5 = np.where(~np.isnan(pndn50_var.flatten()))[0].shape[0]
flt_p9 = np.where(~np.isnan(pndn90_var.flatten()))[0].shape[0]
flt_f5 = np.where(~np.isnan(fndn50_var.flatten()))[0].shape[0]
flt_f9 = np.where(~np.isnan(fndn90_var.flatten()))[0].shape[0]

# make same size as bins
n_fl = np.append(n_fl,0)
n_cl = np.append(n_cl,0)
n_po = np.append(n_po,0)
n_fo = np.append(n_fo,0)
n_p5 = np.append(n_p5,0)
n_p9 = np.append(n_p9,0)
n_f5 = np.append(n_f5,0)
n_f9 = np.append(n_f9,0)

# plot
#plt.ion()

figw = 12
figh = 8
axisfont = 16


fig,ax = plt.subplots(1,1,figsize=[figw,figh])

# plot all bins
ax.plot(bins_fl,n_fl/flt_fl,color='blue',linestyle=':',linewidth=1.5,label='Current day')
ax.plot(bins_cl,n_cl/flt_cl,color='green',linewidth=1.5,label='Ocean only')
ax.plot(bins_po,n_po/flt_po,color='orange',linestyle='-',linewidth=1.5,label='50% N Red.')
ax.plot(bins_p5,n_p5/flt_p5,color='orange',linestyle='--',linewidth=1.5,label='50% N Red. 50% Recy.')
ax.plot(bins_p9,n_p9/flt_p9,color='orange',linestyle=':',linewidth=1.5,label='50% N Red. 90% Recy.')
ax.plot(bins_fo,n_fo/flt_fo,color='gray',linestyle='-',linewidth=1.5,label='85% N Red.')
ax.plot(bins_f5,n_f5/flt_f5,color='gray',linestyle='--',linewidth=1.5,label='85% N Red. 50% Recy.')
ax.plot(bins_f9,n_f9/flt_f9,color='gray',linestyle=':',linewidth=1.5,label='85% N Red. 90% Recy.')

ax.set_yscale('log')
ax.set_xlim([0,45])
ax.set_ylim([1E-6,1E-1])

ax.legend(fontsize=axisfont)
ax.set_xlabel('Surface Chl-a mg m$^{-3}$',fontsize=axisfont)
ax.set_ylabel('PDF',fontsize=axisfont)
ax.set_title(regtitle+' Surface Chl-a PDF',fontsize=axisfont)
ax.tick_params(axis='both',which='major',labelsize=axisfont)
plt.savefig(savepath+savename+'.png',bbox_inches='tight')
