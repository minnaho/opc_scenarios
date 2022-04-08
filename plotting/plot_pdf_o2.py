import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import l2grid
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

plt.ion()

roms_path = '/data/project6/minnaho/opc_scenarios/ext_depth_200/'
savepath = './figs/pdf/'

region_name = 'grid'

varnc = 'var'
varstr = 'O2'

start_year = 1998
end_year = 1998

start_month = 1
end_month = 12

#exp = ['cntrl','l1617','PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
#title_exp = ['CTRL','Loads 16-17','PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']
exp = ['l1617','PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
title_exp = ['Loads 16-17','PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']

if 'cntrl' in exp:
    savename = 'pdf_0_200_'+varstr+'_Y'+str(end_year)+'_cntrl'
else:
    savename = 'pdf_0_200_'+varstr+'_Y'+str(end_year)

filest = 'ext_0_200_'+varstr+'_'

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
mask_nc = l2grid.mask_nc # full grid

mask_ssd[mask_ssd==0] = np.nan
mask_nsd[mask_nsd==0] = np.nan
mask_oc[mask_oc==0] = np.nan
mask_sp[mask_sp==0] = np.nan
mask_sm[mask_sm==0] = np.nan
mask_v[mask_v==0] = np.nan
mask_sb[mask_sb==0] = np.nan
mask_nc[mask_nc==0] = np.nan

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

bin_p = np.load(roms_path+'bin_p_'+str(end_year)+'_'+exp[0]+'.npy')

flt = np.ones((len(exp)))*0
n_p = np.ones((len(exp),bin_p.shape[0]))*0
for e_i in range(len(exp)):
    flt[e_i] = np.load(roms_path+'flt_nonan_'+str(end_year)+'_'+exp[e_i]+'.npy')
    n_p[e_i] = np.load(roms_path+'n_p_count_'+str(end_year)+'_'+exp[e_i]+'.npy')


figw = 12
figh = 8
axisfont = 16
if 'cntrl' in exp:
    cplt = ['green','blue','orange','orange','orange','gray','gray','gray']
    lsty = ['-','-','-','--',':','-','--',':']
else:
    cplt = ['blue','orange','orange','orange','gray','gray','gray']
    lsty = ['-','-','--',':','-','--',':']


fig,ax = plt.subplots(1,1,figsize=[figw,figh])
for e_i in range(len(exp)):
    ax.plot(bin_p,n_p[e_i]/flt[e_i],color=cplt[e_i],linestyle=lsty[e_i],linewidth=1.5,label=title_exp[e_i])

ax.legend(fontsize=axisfont)
ax.set_xlabel('O2 mmol/m3',fontsize=axisfont)
ax.set_ylabel('PDF',fontsize=axisfont)
ax.set_title(regtitle+' O2 PDF '+str(start_year)+'/'+'%02d'%start_month+'-'+str(end_year)+'/'+'%02d'%end_month,fontsize=axisfont)

ax.tick_params(axis='both',which='major',labelsize=axisfont)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.set_xlim([90,300])

plt.savefig(savepath+savename+'.png',bbox_inches='tight')
