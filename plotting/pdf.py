import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import l2grid
import numpy as np
import scipy.stats as stats
from netCDF4 import Dataset,num2date,date2num
import matplotlib.pyplot as plt
from collections import Counter

roms_path = '/data/project6/minnaho/opc_scenarios/ts_int_sli/'
savepath = './figs/'

varnc = 'var'

varname = 'int_100m_biomass'
exp = 'l1617'

# nsd, ssd, oc, sp, sm, v, sb, or none
region_mask = 'sm'

savename = 'pdf_'+varname+'_'+exp

# grid variables
mask_nc = l2grid.mask_nc

# read variables
l1617_nc = Dataset(roms_path+'int_100m_l1617_biomass.nc','r')
cntrl_nc = Dataset(roms_path+'int_100m_cntrl_biomass.nc','r')
fulll_nc = Dataset(roms_path+'int_100m_fulll_biomass.nc','r')

l1617_var = np.array(l1617_nc.variables[varnc])
cntrl_var = np.array(cntrl_nc.variables[varnc])
fulll_var = np.array(fulll_nc.variables[varnc])

l1617_var[l1617_var<0] = np.nan
cntrl_var[cntrl_var<0] = np.nan
fulll_var[fulll_var<0] = np.nan

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

mask_ssd[mask_ssd==0] = np.nan
mask_nsd[mask_nsd==0] = np.nan
mask_oc[mask_oc==0] = np.nan
mask_sp[mask_sp==0] = np.nan
mask_sm[mask_sm==0] = np.nan
mask_v[mask_v==0] = np.nan
mask_sb[mask_sb==0] = np.nan

if region_mask == 'ssd':
    mask_mult = mask_ssd
    regtitle = 'South San Diego'
if region_mask == 'nsd':
    mask_mult = mask_nsd
    regtitle = 'North San Diego'
if region_mask == 'oc':
    mask_mult = mask_oc
    regtitle = 'Orange County'
if region_mask == 'sp':
    mask_mult = mask_sp
    regtitle = 'San Pedro'
if region_mask == 'sm':
    mask_mult = mask_sm
    regtitle = 'Santa Monica Bay'
if region_mask == 'v':
    mask_mult = mask_v
    regtitle = 'Ventura'
if region_mask == 'sb':
    mask_mult = mask_sb
    regtitle = 'Santa Barbara'
if region_mask == 'none':
    mask_mult = mask_nc
    regtitle = 'SCB'

# multiply by region mask
# comment out for full region
#l1617_var = mask_mult*l1617_var
#cntrl_var = mask_mult*cntrl_var
#fulll_var = mask_mult*fulll_var

# stats
nbins = 500
n_c,bins_c,patch_c = plt.hist(cntrl_var.flatten(),nbins)
n_f,bins_f,patch_f = plt.hist(fulll_var.flatten(),nbins)
n_l,bins_l,patch_l = plt.hist(l1617_var.flatten(),nbins)

# non nan values
flt_c = np.where(~np.isnan(cntrl_var.flatten()))[0].shape[0]
flt_f = np.where(~np.isnan(fulll_var.flatten()))[0].shape[0]
flt_l = np.where(~np.isnan(l1617_var.flatten()))[0].shape[0]
#flt = cntrl_var.flatten().shape[0]

#l1617_stat = stats.gaussian_kde(l1617_var)
#cntrl_stat = stats.gaussian_kde(cntrl_var)
#fulll_stat = stats.gaussian_kde(fulll_var)

# make same size as bins
n_l = np.append(n_l,0)
n_c = np.append(n_c,0)
n_f = np.append(n_f,0)

# plot
plt.ion()

figw = 12
figh = 8
axisfont = 16


fig,ax = plt.subplots(1,1,figsize=[figw,figh])

# remove first bin that has 0 values 
# divide by shape to get probability
#ax.plot(bins_c[1:],n_c[1:]/flt_c,color='green',linewidth=1.5,label='CTRL')
#ax.plot(bins_f[1:],n_f[1:]/flt_f,color='black',linewidth=1.5,label='FULL')
#ax.plot(bins_l[1:],n_l[1:]/flt_l,color='blue',linewidth=1.5,label='16-17 loads')

ax.plot(bins_c,n_c/flt_c,color='green',linewidth=1.5,label='CTRL')
ax.plot(bins_f,n_f/flt_f,color='black',linewidth=1.5,label='FULL')
ax.plot(bins_l,n_l/flt_l,color='blue',linewidth=1.5,label='16-17 loads')

ax.set_xlim([-10,1000])

ax.legend(fontsize=axisfont)
ax.set_xlabel('integrated biomass 100 m mmol/m2',fontsize=axisfont)
ax.set_ylabel('PDF',fontsize=axisfont)
ax.set_title(regtitle+' integrated biomass PDF',fontsize=axisfont)
ax.tick_params(axis='both',which='major',labelsize=axisfont)
plt.savefig(savepath+savename,bbox_inches='tight')
