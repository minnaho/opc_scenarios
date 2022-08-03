# bar plot of difference only
import sys
import os
sys.path.append('/data/project3/minnaho/global/')
import l2grid
import numpy as np
from netCDF4 import Dataset,num2date
import glob as glob
import matplotlib.pyplot as plt

plt.ion()

savepath = './figs/scatter/'
region_name = 'offshore'

# ROMS output location
outpath = '/data/project6/minnaho/opc_scenarios/bgc_flux/'

# roms var
var_name = 'npp' 
var_nc = 'var_int' 
cblabel = 'mmol m$^{-2}$ d$^{-1}$'

filest = 'int_avg_100m_50m_'

#year_month = 'Y1998_M04_06'
year_month1 = 'fullts'
year_month2 = 'fullts'

# scenario names 
exp = ['PNDN_only',
       'FNDN_only',
       'pndn50', 
       'pndn90', 
       'fndn50',
       'fndn90']

title_exp = ['50% N Reduction',
             '85% N Reduction',
             '50% N Reduction\n50% Recycle',
             '50% N Reduction\n90% Recycle',
             '85% N Reduction\n50% Recycle',
             '85% N Reduction\n90% Recycle']

# mask
mask_nc = l2grid.mask_nc

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
mask_cst[mask_cst==0] = np.nan

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
    mask_temp = np.copy(mask_nc)
    mask_temp[:,:20] = np.nan
    mask_temp[:20,:] = np.nan
    mask_temp[-20:,:] = np.nan
    mask_mult = mask_temp
    regtitle = 'full SCB'
if region_name == 'offshore':
    # do opposite of coastal band mask
    mask_temp = np.copy(mask_cst)
    mask_temp[mask_temp==1] = 2
    mask_temp[np.isnan(mask_temp)] = 1
    mask_temp[mask_temp==2] = 0
    mask_temp = mask_temp*mask_nc
    mask_temp[:,:20] = np.nan
    mask_temp[:20,:] = np.nan
    mask_temp[-20:,:] = np.nan
    mask_mult = mask_temp
    regtitle = 'offshore'

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

# seconds to days
s2d = 86400

filest1 = 'avg_'+year_month1+'_int_avg_100m_50m_'
filest2 = 'avg_'+year_month2+'_int_avg_100m_50m_'
fileen = '_'+var_name+'.nc'


# get cntrl and loads1617 to compare to
cntrl_var_old = np.squeeze(Dataset(outpath+'avg_'+year_month1+'_int_avg_100m_50m_cntrl_initap'+fileen,'r').variables[var_nc])*mask_mult
cntrl_var_old[cntrl_var_old==0] = np.nan
cntrl_mean_old = np.nanmean(cntrl_var_old)*s2d

cntrl_var_new = np.squeeze(Dataset(outpath+'avg_'+year_month2+'_int_avg_100m_50m_cntrl_2012_2017'+fileen,'r').variables[var_nc])*mask_mult
cntrl_var_new[cntrl_var_new==0] = np.nan
cntrl_mean_new = np.nanmean(cntrl_var_new)*s2d


fulll_var_old = np.squeeze(Dataset(outpath+'avg_'+year_month1+'_int_avg_100m_50m_loads1617'+fileen,'r').variables[var_nc])*mask_mult
fulll_var_old[fulll_var_old==0] = np.nan
fulll_mean_old = np.nanmean(fulll_var_old)*s2d

fulll_var_new = np.squeeze(Dataset(outpath+'avg_'+year_month2+'_int_avg_100m_50m_fulll_2012_2017'+fileen,'r').variables[var_nc])*mask_mult
fulll_var_new[fulll_var_new==0] = np.nan
fulll_mean_new = np.nanmean(fulll_var_new)*s2d

clim = np.ones((len(exp),12))*np.nan
cstd = np.ones((len(exp),12))*np.nan
# per month
for m_i in range(1,13):
    for e_i in range(len(exp)):
        flist = glob.glob(outpath+filest+exp[e_i]+'*'+var_name+'*M%02d'%m_i+'.nc')
        temparr = np.ones((len(flist),mask_nc.shape[0],mask_nc.shape[1]))*np.nan
        for f_i in range(len(flist)):
            temprd = np.squeeze(Dataset(flist[f_i],'r').variables[var_nc])*mask_mult
            temprd[temprd==0] = np.nan
            if 'real' in flist[f_i]:
                cntrlm = np.squeeze(Dataset(outpath+filest+'cntrl_2012_2017_'+var_name+'_'+flist[f_i][flist[f_i].index('Y'):],'r').variables[var_nc])*mask_mult
            else:
                cntrlm = np.squeeze(Dataset(outpath+filest+'cntrl_initap_'+var_name+'_'+flist[f_i][flist[f_i].index('Y'):],'r').variables[var_nc])*mask_mult
            cntrlm[cntrlm==0] = np.nan
            temparr[f_i] = temprd - cntrlm
        cstd[e_i,m_i-1] = np.nanstd(np.nanmean(temparr,axis=(1,2)),axis=0)
        clim[e_i,m_i-1] = np.nanmean(temparr)

clim = clim*s2d
cstd = cstd*s2d

# full ANTH climatology
clim_fulll = np.ones((12))*np.nan
cstd_fulll = np.ones((12))*np.nan
# per month
for m_i in range(1,13):
    flist1 = glob.glob(outpath+filest+'fulll*'+var_name+'*M%02d'%m_i+'.nc')
    flist2 = glob.glob(outpath+filest+'loads1617*'+var_name+'*M%02d'%m_i+'.nc')
    flist_fulll = flist1+flist2
    temparr_fulll = np.ones((len(flist_fulll),mask_nc.shape[0],mask_nc.shape[1]))*np.nan
    
    for f_i in range(len(flist_fulll)):
        temprd_fulll = np.squeeze(Dataset(flist_fulll[f_i],'r').variables[var_nc])*mask_mult
        temprd_fulll[temprd_fulll==0] = np.nan
        if '2012_2017' in flist_fulll[f_i]:
            cntrlm = np.squeeze(Dataset(outpath+filest+'cntrl_2012_2017_'+var_name+'_'+flist_fulll[f_i][flist_fulll[f_i].index('Y'):],'r').variables[var_nc])*mask_mult
        else:
            cntrlm = np.squeeze(Dataset(outpath+filest+'cntrl_initap_'+var_name+'_'+flist_fulll[f_i][flist_fulll[f_i].index('Y'):],'r').variables[var_nc])*mask_mult
        cntrlm[cntrlm==0] = np.nan
        temparr_fulll[f_i] = temprd_fulll - cntrlm

    cstd_fulll[m_i-1] = np.nanstd(np.nanmean(temparr_fulll,axis=(1,2)),axis=0)
    clim_fulll[m_i-1] = np.nanmean(temparr_fulll)

clim_fulll = clim_fulll*s2d
cstd_fulll = cstd_fulll*s2d


figw = 16
figh = 4

months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

axis_tick_size = 14


for n_i in range(len(exp)):
    fig,ax = plt.subplots(1,1,figsize=[figw,figh])
    ax.bar(range(1,13),clim[n_i],yerr=cstd[n_i],color='white',edgecolor='k',capsize=4)

    if region_name == 'coast':
        ax.set_ylim(bottom=-5,top=30)
    if region_name == 'grid':
        ax.set_ylim(bottom=-10,top=10)
    if region_name == 'offshore':
        ax.set_ylim(bottom=-8,top=8)
    #ax.set_ylim(bottom=0)
    ax.set_xticks(range(1,13))
    ax.plot(range(1,13),clim_fulll)
    ax.fill_between(range(1,13),clim_fulll+cstd_fulll,clim_fulll-cstd_fulll,alpha=0.3)
    
    ax.set_xticklabels(months,fontsize=12)
    ax.set_ylabel('Integrated NPP 100 m '+cblabel,fontsize=axis_tick_size)
    ax.set_title(title_exp[n_i],fontsize=12)
    ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)
    
    savename = var_name+'_'+exp[n_i]+'_clim_4years_diff'+region_name+'.png'
    fig.savefig(savepath+savename,bbox_inches='tight')
    print(savename)
    #plt.close()

