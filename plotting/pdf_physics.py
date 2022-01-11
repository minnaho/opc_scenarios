import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import l2grid
import numpy as np
import scipy.stats as stats
from netCDF4 import Dataset,num2date,date2num
import matplotlib.pyplot as plt
from collections import Counter

roms_path = '/data/project6/minnaho/opc_scenarios/physics/'
savepath = './figs/pdf/'

varnc = 'u'

varname = 'u'
exp = 'Y1998M02'

# nsd, ssd, oc, sp, sm, v, sb, coast, or scb 
region_name = 'scb'

savename = 'pdf_'+varname+'_'+exp+'_'+region_name

# grid variables
mask_nc = l2grid.mask_nc

dtstr = 'Y1998M02'

# read variables
l1617_nc = Dataset(roms_path+'l2_scb_avg.'+dtstr+'_uv_concat_loads1617.nc','r')
cntrl_nc = Dataset(roms_path+'l2_scb_avg.'+dtstr+'_uv_concat_pndn90.nc','r')
fulll_nc = Dataset(roms_path+'l2_scb_avg.'+dtstr+'_uv_concat_fndn90.nc','r')
exp01_nc = Dataset(roms_path+'l2_scb_avg.'+dtstr+'_uv_concat_PNDN_only.nc','r')
exp02_nc = Dataset(roms_path+'l2_scb_avg.'+dtstr+'_uv_concat_FNDN_only.nc','r')
exp03_nc = Dataset(roms_path+'l2_scb_avg.'+dtstr+'_uv_concat_pndn50.nc','r')
#exp04_nc = Dataset(roms_path+'int_100m_pndn90_biomass_'+dtstr+'.nc','r')

l1617_var = np.array(l1617_nc.variables[varnc])
cntrl_var = np.array(cntrl_nc.variables[varnc])
fulll_var = np.array(fulll_nc.variables[varnc])
exp01_var = np.array(exp01_nc.variables[varnc])
exp02_var = np.array(exp02_nc.variables[varnc])
exp03_var = np.array(exp03_nc.variables[varnc])



# remove all nan values
l1617_var[l1617_var>1E10] = np.nan
cntrl_var[cntrl_var>1E10] = np.nan
fulll_var[fulll_var>1E10] = np.nan
exp01_var[exp01_var>1E10] = np.nan
exp02_var[exp02_var>1E10] = np.nan
exp03_var[exp03_var>1E10] = np.nan

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
if region_name == 'scb': # full L2 grid
    mask_mult = mask_nc
    regtitle = 'SCB'
if region_name == 'coast':
    mask_mult = mask_cst
    regtitle = '15 km coast'

#l1617_var = np.ones((l1617_var_nomask.shape[0],l1617_var_nomask.shape[1],
#                     l1617_var_nomask.shape[2],l1617_var_nomask.shape[3]))*np.nan
#
#for t_i in range(l1617_var.shape[0]):
#    l1617_var[t_i] = l1617_var_nomask[t_i]*mask_mult
#    #cntrl_var[t_i] = cntrl_var[t_i]*mask_mult
#    #fulll_var[t_i] = fulll_var[t_i]*mask_mult
#    #exp01_var[t_i] = exp01_var[t_i]*mask_mult
#    #exp02_var[t_i] = exp02_var[t_i]*mask_mult
#    #exp03_var[t_i] = exp03_var[t_i]*mask_mult


# stats
nbins = 500
n_c,bins_c,patch_c = plt.hist(cntrl_var.flatten(),nbins)
n_f,bins_f,patch_f = plt.hist(fulll_var.flatten(),nbins)
n_l,bins_l,patch_l = plt.hist(l1617_var.flatten(),nbins)
n_1,bins_1,patch_1 = plt.hist(exp01_var.flatten(),nbins)
n_2,bins_2,patch_2 = plt.hist(exp02_var.flatten(),nbins)
n_3,bins_3,patch_3 = plt.hist(exp03_var.flatten(),nbins)

# non nan values
flt_c = np.where(~np.isnan(cntrl_var.flatten()))[0].shape[0]
flt_f = np.where(~np.isnan(fulll_var.flatten()))[0].shape[0]
flt_l = np.where(~np.isnan(l1617_var.flatten()))[0].shape[0]
flt_1 = np.where(~np.isnan(exp01_var.flatten()))[0].shape[0]
flt_2 = np.where(~np.isnan(exp02_var.flatten()))[0].shape[0]
flt_3 = np.where(~np.isnan(exp03_var.flatten()))[0].shape[0]
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

# plot
plt.ion()

figw = 12
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
ax.plot(bins_l,n_l/flt_l,color='blue',linestyle=':',linewidth=1.5,label='16-17 loads')
ax.plot(bins_1,n_1/flt_1,color='orange',linestyle='--',linewidth=1.5,label='PNDN only')
ax.plot(bins_3,n_3/flt_3,color='purple',linestyle='dashdot',linewidth=1.5,label='PNDN 50')
ax.plot(bins_c,n_c/flt_c,color='green',linewidth=1.5,label='PNDN 90')
ax.plot(bins_2,n_2/flt_2,color='gray',linestyle='--',linewidth=1.5,label='FNDN only')
ax.plot(bins_f,n_f/flt_f,color='r',linestyle=':',linewidth=1.5,label='FNDN 90')

#ax.set_xlim([-10,500])

ax.legend(fontsize=axisfont)
ax.set_xlabel(varname+' m/s',fontsize=axisfont)
ax.set_ylabel('PDF',fontsize=axisfont)
ax.set_title(regtitle+' '+varname+' PDF',fontsize=axisfont)
ax.tick_params(axis='both',which='major',labelsize=axisfont)
plt.savefig(savepath+savename+'.png',bbox_inches='tight')

print('skew 16-17',stats.skew(n_l/flt_l))
print('skew PNDN',stats.skew(n_1/flt_1))
print('skew PNDN50',stats.skew(n_3/flt_3))
print('skew PNDN90',stats.skew(n_c/flt_c))
print('skew FNDN',stats.skew(n_2/flt_2))
print('skew FNDN90',stats.skew(n_f/flt_f))

print('variance 16-17',np.var(n_l/flt_l))
print('variance PNDN',np.var(n_1/flt_1))
print('variance PNDN50',np.var(n_3/flt_3))
print('variance PNDN90',np.var(n_c/flt_c))
print('variance FNDN',np.var(n_2/flt_2))
print('variance FNDN90',np.var(n_f/flt_f))

