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
savepath = './figs/pdf/'

varnc = 'var'

varname = 'int_100m_biomass'
exp = '1998'

# nsd, ssd, oc, sp, sm, v, sb, or grid or coast (15 km)
region_name = 'grid'

# percentile of cntrl to take to compare frequencies
perc = 75

#savename = 'pdf_fullts_'+varname+'_0bin_'+region_name
savename = 'freq_pdf_'+varname+'_'+exp+'_0bin_'+region_name+'_'+str(perc)

# grid variables
mask_nc = l2grid.mask_nc

dtstr = 'Y1998M01_12'
#dtstr = 'Y1999M01_11'


# read variables
l1617_nc = Dataset(roms_path+'concat_int_100m_l1617_biomass_'+dtstr+'.nc','r')
cntrl_nc = Dataset(roms_path+'concat_int_100m_cntrl_biomass_'+dtstr+'.nc','r')
fulll_nc = Dataset(roms_path+'concat_int_100m_fndn90_biomass_'+dtstr+'.nc','r')
exp01_nc = Dataset(roms_path+'concat_int_100m_PNDN_only_biomass_'+dtstr+'.nc','r')
exp02_nc = Dataset(roms_path+'concat_int_100m_FNDN_only_biomass_'+dtstr+'.nc','r')
exp03_nc = Dataset(roms_path+'concat_int_100m_pndn50_biomass_'+dtstr+'.nc','r')
exp04_nc = Dataset(roms_path+'concat_int_100m_pndn90_biomass_'+dtstr+'.nc','r')
exp05_nc = Dataset(roms_path+'concat_int_100m_fndn50_biomass_'+dtstr+'.nc','r')


l1617_var = np.array(l1617_nc.variables[varnc])
cntrl_var = np.array(cntrl_nc.variables[varnc])
fulll_var = np.array(fulll_nc.variables[varnc])
exp01_var = np.array(exp01_nc.variables[varnc])
exp02_var = np.array(exp02_nc.variables[varnc])
exp03_var = np.array(exp03_nc.variables[varnc])
exp04_var = np.array(exp04_nc.variables[varnc])
exp05_var = np.array(exp05_nc.variables[varnc])

# remove all negative and 0 values
l1617_var[l1617_var<=0] = np.nan
cntrl_var[cntrl_var<=0] = np.nan
fulll_var[fulll_var<=0] = np.nan
exp01_var[exp01_var<=0] = np.nan
exp02_var[exp02_var<=0] = np.nan
exp03_var[exp03_var<=0] = np.nan
exp04_var[exp04_var<=0] = np.nan
exp05_var[exp05_var<=0] = np.nan

# remove all negative values
#l1617_var[l1617_var<0] = np.nan
#cntrl_var[cntrl_var<0] = np.nan
#fulll_var[fulll_var<0] = np.nan
#exp01_var[exp01_var<0] = np.nan
#exp02_var[exp02_var<0] = np.nan
#exp03_var[exp03_var<0] = np.nan

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
    regtitle = '15 km coast'

l1617_var = l1617_var*mask_mult
cntrl_var = cntrl_var*mask_mult
fulll_var = fulll_var*mask_mult
exp01_var = exp01_var*mask_mult
exp02_var = exp02_var*mask_mult
exp03_var = exp03_var*mask_mult
exp04_var = exp04_var*mask_mult
exp05_var = exp05_var*mask_mult

# stats
nbins = 500
n_c,bins_c,patch_c = plt.hist(cntrl_var.flatten(),nbins)
n_f,bins_f,patch_f = plt.hist(fulll_var.flatten(),nbins)
n_l,bins_l,patch_l = plt.hist(l1617_var.flatten(),nbins)
n_1,bins_1,patch_1 = plt.hist(exp01_var.flatten(),nbins)
n_2,bins_2,patch_2 = plt.hist(exp02_var.flatten(),nbins)
n_3,bins_3,patch_3 = plt.hist(exp03_var.flatten(),nbins)
n_4,bins_4,patch_4 = plt.hist(exp04_var.flatten(),nbins)
n_5,bins_5,patch_5 = plt.hist(exp05_var.flatten(),nbins)

# non nan values
flt_c = np.where(~np.isnan(cntrl_var.flatten()))[0].shape[0]
flt_f = np.where(~np.isnan(fulll_var.flatten()))[0].shape[0]
flt_l = np.where(~np.isnan(l1617_var.flatten()))[0].shape[0]
flt_1 = np.where(~np.isnan(exp01_var.flatten()))[0].shape[0]
flt_2 = np.where(~np.isnan(exp02_var.flatten()))[0].shape[0]
flt_3 = np.where(~np.isnan(exp03_var.flatten()))[0].shape[0]
flt_4 = np.where(~np.isnan(exp04_var.flatten()))[0].shape[0]
flt_5 = np.where(~np.isnan(exp05_var.flatten()))[0].shape[0]
#flt = cntrl_var.flatten().shape[0]

#l1617_stat = stats.gaussian_kde(l1617_var)
#cntrl_stat = stats.gaussian_kde(cntrl_var)
#fulll_stat = stats.gaussian_kde(fulll_var)

# make same size as bins
n_l = np.append(n_l,0)
n_c = np.append(n_c,0)
n_f = np.append(n_f,0)
n_1 = np.append(n_1,0)
n_2 = np.append(n_2,0)
n_3 = np.append(n_3,0)
n_4 = np.append(n_4,0)
n_5 = np.append(n_5,0)

# calculate frequency above xth percentile of cntrl

cntrlflat = cntrl_var.flatten()
cntrlperc = np.nanpercentile(cntrlflat,perc)

cntrl_ind = np.where(bins_c>cntrlperc)
l1617_ind = np.where(bins_l>cntrlperc)
pndnon_ind = np.where(bins_1>cntrlperc)
pndn50_ind = np.where(bins_3>cntrlperc)
pndn90_ind = np.where(bins_4>cntrlperc)
fndnon_ind = np.where(bins_2>cntrlperc)
fndn50_ind = np.where(bins_5>cntrlperc)
fndn90_ind = np.where(bins_f>cntrlperc)

cntrl_freq =  np.nansum((n_c/flt_c)[cntrl_ind])
l1617_freq =  np.nansum((n_l/flt_c)[l1617_ind])
pndnon_freq = np.nansum((n_1/flt_c)[pndnon_ind])
pndn50_freq = np.nansum((n_3/flt_c)[pndn50_ind])
pndn90_freq = np.nansum((n_4/flt_c)[pndn90_ind])
fndnon_freq = np.nansum((n_2/flt_c)[fndnon_ind])
fndn50_freq = np.nansum((n_5/flt_c)[fndn50_ind])
fndn90_freq = np.nansum((n_f/flt_c)[fndn90_ind])

# plot
#plt.ion()

figw = 14
figh = 8
axisfont = 16


fig,ax = plt.subplots(1,1,figsize=[figw,figh])

# remove first bin that has 0 values 
# divide by shape to get probability
#ax.plot(bins_c[1:],n_c[1:]/flt_c,color='green',linewidth=1.5,label='CTRL')
#ax.plot(bins_l[1:],n_l[1:]/flt_l,color='blue',linestyle=':',linewidth=1.5,label='16-17 loads')
#ax.plot(bins_1[1:],n_1[1:]/flt_1,color='orange',linestyle='--',linewidth=1.5,label='PNDN only')
#ax.plot(bins_3[1:],n_3[1:]/flt_3,color='lightblue',linestyle='--',linewidth=1.5,label='PNDN 50')
#ax.plot(bins_2[1:],n_2[1:]/flt_2,color='gray',linestyle='--',linewidth=1.5,label='FNDN only')
#ax.plot(bins_f[1:],n_f[1:]/flt_f,color='black',linewidth=1.5,label='FNDN 90')

# plot all bins
title_exp = ['CTRL','Loads 16-17','PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']

n_i = range(len(title_exp))

ax.bar(n_i[0],cntrl_freq)
ax.bar(n_i[1],l1617_freq)
ax.bar(n_i[2],pndnon_freq)
ax.bar(n_i[3],pndn50_freq)
ax.bar(n_i[4],pndn90_freq)
ax.bar(n_i[5],fndnon_freq)
ax.bar(n_i[6],fndn50_freq)
ax.bar(n_i[7],fndn90_freq)

ax.set_xticks(range(len(title_exp)))
ax.set_xticklabels(title_exp)

ax.set_ylabel('Integrated Biomass Frequency > '+str(perc)+'th CTRL percentile',fontsize=axisfont)
#ax.set_title(regtitle+' integrated biomass PDF',fontsize=axisfont)
ax.tick_params(axis='both',which='major',labelsize=axisfont)
plt.savefig(savepath+savename+'.png',bbox_inches='tight')
