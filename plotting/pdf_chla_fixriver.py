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

savename = 'pdf_'+varname+'_'+region_name+'_log_fixriver'

# grid variables
mask_nc = l2grid.mask_nc


# read variables
fulll_nc = Dataset(roms_path+'concat_totchla_fulll_2012_2017.nc','r')
cntrl_nc = Dataset(roms_path+'concat_totchla_cntrl_initap_realistic.nc','r')
pndnon_nc = Dataset(roms_path+'concat_totchla_PNDN_only_realistic.nc','r')
fndnon_nc = Dataset(roms_path+'concat_totchla_FNDN_only_realistic.nc','r')
pndn50_nc = Dataset(roms_path+'concat_totchla_pndn50_fixriver.nc','r')
pndn90_nc = Dataset(roms_path+'concat_totchla_pndn90_fixriver.nc','r')
fndn50_nc = Dataset(roms_path+'concat_totchla_fndn50_fixriver.nc','r')
fndn90_nc = Dataset(roms_path+'concat_totchla_fndn90_fixriver.nc','r')

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

fulll_mask = np.squeeze(fulll_nc.variables[varnc])*mask_mult
cntrl_mask = np.squeeze(cntrl_nc.variables[varnc])*mask_mult
pndnon_mask = np.squeeze(pndnon_nc.variables[varnc])*mask_mult
fndnon_mask = np.squeeze(fndnon_nc.variables[varnc])*mask_mult
pndn50_mask = np.squeeze(pndn50_nc.variables[varnc])*mask_mult
pndn90_mask = np.squeeze(pndn90_nc.variables[varnc])*mask_mult
fndn50_mask = np.squeeze(fndn50_nc.variables[varnc])*mask_mult
fndn90_mask = np.squeeze(fndn90_nc.variables[varnc])*mask_mult

fulll_mask = fulll_mask[fulll_mask>0]
cntrl_mask = cntrl_mask[cntrl_mask>0]
pndnon_mask = pndnon_mask[pndnon_mask>0]
fndnon_mask = fndnon_mask[fndnon_mask>0]
pndn50_mask = pndn50_mask[pndn50_mask>0]
pndn90_mask = pndn90_mask[pndn90_mask>0]
fndn50_mask = fndn50_mask[fndn50_mask>0]
fndn90_mask = fndn90_mask[fndn90_mask>0]

fulll_var = fulll_mask
cntrl_var = cntrl_mask
pndnon_var =pndnon_mask
fndnon_var =fndnon_mask
pndn50_var =pndn50_mask
pndn90_var =pndn90_mask
fndn50_var =fndn50_mask
fndn90_var =fndn90_mask

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
plt.ion()

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
ax.set_xlim([0,20])
ax.set_ylim([1E-6,1E-1])

ax.legend(fontsize=axisfont)
ax.set_xlabel('Surface Chl-a mg m$^{-3}$',fontsize=axisfont)
ax.set_ylabel('PDF',fontsize=axisfont)
ax.set_title(regtitle+' Surface Chl-a PDF',fontsize=axisfont)
ax.tick_params(axis='both',which='major',labelsize=axisfont)
plt.savefig(savepath+savename+'.png',bbox_inches='tight')
