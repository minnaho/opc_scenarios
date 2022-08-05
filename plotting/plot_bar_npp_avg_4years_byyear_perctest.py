# bar plot of difference only
# plot of each year: 1998,1999,2016,2017
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

# months to average over for offshore/full bight
# 1997:1997M11-1998M10,1998M11-1999M10,2015M11-2016M10,2016M11-2017M10
# fulltime:avg per year of 1998,1999,2016,2017
# junoct - jun to oct
# springsummer = apr-sep 

# coast time period
ctimep = '1997'
# grid/offshore time period
gtimep = '1997'


# ROMS output location
outpath = '/data/project6/minnaho/opc_scenarios/bgc_flux/'

# roms var
var_name = 'npp' 
var_nc = 'var_int' 
cblabel = 'mmol m$^{-2}$ d$^{-1}$'

filest = 'int_avg_100m_50m_'

yearlist1 = [1998,1999]
yearlist2 = [2016,2017]

# scenario names 
exp1 = ['PNDN_only',
        'pndn50',
        'pndn90',
       'FNDN_only',
       'fndn50',
       'fndn90'
       ]

exp2 = ['PNDN_only_realistic',
        'pndn50_realistic',
        'pndn90_realistic',
        'FNDN_only_realistic'
       ]

title_exp1 = ['50% N Red.',
             '50% N Red.\n50% Recy.',
             '50% N Red.\n90% Recy.',
             '85% N Red.',
             '85% N Red.\n50% Recy.',
             '85% N Red.\n90% Recy.'
             ]

title_exp2 = ['50% N Red.',
              '50% N Red.\n50% Recy.',
              '50% N Red.\n90% Recy.',
              '85% N Red.'
             ]

#title_exp = ['50% N Reduction\n98-99',
#             '16-17',
#             '50% N Reduction\n50% Recycle\n98-99',
#             '16-17',
#             '50% N Reduction\n90% Recycle\n98-99',
#             '16-17',
#             '85% N Reduction\n98-99',
#             '16-17',
#             '85% N Reduction\n50% Recycle',
#             '85% N Reduction\n90% Recycle']

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
    regtitle = '15 km Coast'
if region_name == 'grid':
    mask_temp = np.copy(mask_nc)
    mask_temp[:,:20] = np.nan
    mask_temp[:20,:] = np.nan
    mask_temp[-20:,:] = np.nan
    mask_mult = mask_temp
    regtitle = 'Bightwide'

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

if region_name == 'offshore':
    # do opposite of coastal band mask
    '''
    mask_temp = np.copy(mask_cst)
    mask_temp[mask_temp==1] = 2
    mask_temp[np.isnan(mask_temp)] = 1
    mask_temp[mask_temp==2] = 0
    mask_temp = mask_temp*mask_nc
    mask_temp[:,:20] = np.nan
    mask_temp[:20,:] = np.nan
    mask_temp[-20:,:] = np.nan
    mask_mult = mask_temp
    regtitle = 'Offshore'
    '''
    # faycal's offshore mask7+mask9
    mask7[mask9==1] = 1 
    #mask7[mask8==1] = 1 
    mask_mult = mask7

# seconds to days
s2d = 86400

fileen = '_'+var_name+'.nc'

# 1998-1999
avganth1 = np.ones((len(yearlist1)))*np.nan
avgp1 = np.ones((len(exp1),len(yearlist1)))*np.nan
stdanth1 = np.ones((len(yearlist1)))*np.nan
stdp1 = np.ones((len(exp1),len(yearlist1)))*np.nan
for y_i1 in range(len(yearlist1)):
    for e_i1 in range(len(exp1)):
        # get all months
        if region_name == 'coast':
            if ctimep == 'fulltime':
                flist = glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1])+'M??.nc')
            if ctimep == 'springsummer':
                flist = glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1])+'M0[4-9].nc')

            # try 1997M11-1998M10 and 1997M11-1998M10
            if ctimep == '1997':
                flist = glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1]-1)+'M1[1-2].nc')+glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1])+'M0[1-9].nc')+glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1])+'M10.nc')

        # only get certain months for grid/offshore 
        elif region_name == 'grid' or region_name == 'offshore':
            if gtimep == 'fulltime':
                flist = glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1])+'M??.nc')
            if gtimep == 'springsummer':
                flist = glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1])+'M0[4-9].nc')
            if gtimep == 'junoct':
                flist = glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1])+'M0[6-9].nc')+glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1])+'M10.nc')
            if gtimep == 'sepnov':
                flist = glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1])+'M09.nc')+glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1])+'M1[0-1].nc')

            # try 1997M11-1998M10 and 1997M11-1998M10
            if gtimep == '1997':
                flist = glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1]-1)+'M1[1-2].nc')+glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1])+'M0[1-9].nc')+glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1])+'M10.nc')

        temparr = np.ones((len(flist),mask_nc.shape[0],mask_nc.shape[1]))*np.nan
        for f_i in range(len(flist)):
            temprd = np.squeeze(Dataset(flist[f_i],'r').variables[var_nc])
            cntrlm = np.squeeze(Dataset(outpath+filest+'cntrl_initap_'+var_name+'_'+flist[f_i][flist[f_i].index('Y'):],'r').variables[var_nc])
            temparr[f_i] = (temprd - cntrlm)*mask_mult*s2d
            temparr[temparr==0] = np.nan

        avgp1[e_i1,y_i1] = np.nanmean(temparr) 
        stdp1[e_i1,y_i1] = np.nanstd(np.nanmean(temparr,axis=(1,2)))

    # now calculate for anth run
    if region_name == 'coast':
        if ctimep == 'fulltime':
            flist = glob.glob(outpath+filest+'loads1617_'+var_name+'_Y'+str(yearlist1[y_i1])+'M??.nc')
        if ctimep == 'springsummer':
            flist = glob.glob(outpath+filest+'loads1617_'+var_name+'_Y'+str(yearlist1[y_i1])+'M0[4-9].nc')

        # try 1997M11-1998M10 and 1997M11-1998M10
        if ctimep == '1997':
            flist = glob.glob(outpath+filest+'loads1617_'+var_name+'_Y'+str(yearlist1[y_i1]-1)+'M1[1-2].nc')+glob.glob(outpath+filest+'loads1617_'+var_name+'_Y'+str(yearlist1[y_i1])+'M0[1-9].nc')+glob.glob(outpath+filest+'loads1617_'+var_name+'_Y'+str(yearlist1[y_i1])+'M10.nc')

    # only get certain months for grid/offshore 
    elif region_name == 'grid' or region_name == 'offshore':
        if gtimep == 'fulltime':
            flist = glob.glob(outpath+filest+'loads1617_'+var_name+'_Y'+str(yearlist1[y_i1])+'M??.nc')
        if gtimep == 'springsummer':
            flist = glob.glob(outpath+filest+'loads1617_'+var_name+'_Y'+str(yearlist1[y_i1])+'M0[4-9].nc')
        if gtimep == 'junoct':
            flist = glob.glob(outpath+filest+'loads1617_'+var_name+'_Y'+str(yearlist1[y_i1])+'M0[6-9].nc')+glob.glob(outpath+filest+'loads1617_'+var_name+'_Y'+str(yearlist1[y_i1])+'M10.nc')
        if gtimep == 'sepnov':
            flist = glob.glob(outpath+filest+'loads1617_'+var_name+'_Y'+str(yearlist1[y_i1])+'M09.nc')+glob.glob(outpath+filest+'loads1617_'+var_name+'_Y'+str(yearlist1[y_i1])+'M1[0-1].nc')

        # try 1997M11-1998M10 and 1997M11-1998M10
        if gtimep == '1997':
            flist = glob.glob(outpath+filest+'loads1617_'+var_name+'_Y'+str(yearlist1[y_i1]-1)+'M1[1-2].nc')+glob.glob(outpath+filest+'loads1617_'+var_name+'_Y'+str(yearlist1[y_i1])+'M0[1-9].nc')+glob.glob(outpath+filest+'loads1617_'+'_'+var_name+'_Y'+str(yearlist1[y_i1])+'M10.nc')

    temparr = np.ones((len(flist),mask_nc.shape[0],mask_nc.shape[1]))*np.nan
    for f_i in range(len(flist)):
        temprd = np.squeeze(Dataset(flist[f_i],'r').variables[var_nc])
        cntrlm = np.squeeze(Dataset(outpath+filest+'cntrl_initap_'+var_name+'_'+flist[f_i][flist[f_i].index('Y'):],'r').variables[var_nc])
        temparr[f_i] = (temprd - cntrlm)*mask_mult*s2d
        temparr[temparr==0] = np.nan
    avganth1[y_i1] = np.nanmean(temparr) 
    stdanth1[y_i1] = np.nanstd(np.nanmean(temparr,axis=(1,2)))
    

# 2016-2017
avganth2 = np.ones((len(yearlist2)))*np.nan
avgp2 = np.ones((len(exp2),len(yearlist2)))*np.nan
stdanth2 = np.ones((len(yearlist2)))*np.nan
stdp2 = np.ones((len(exp2),len(yearlist2)))*np.nan
for y_i2 in range(len(yearlist2)):
    for e_i2 in range(len(exp2)):
        if region_name == 'coast':
            if ctimep == 'fulltime':
                flist = glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2])+'M??.nc')
            if ctimep == 'springsummer':
                flist = glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2])+'M0[4-9].nc')

            # try 2015M11-2016M10 and 2016M11-2017M10
            if ctimep == '1997':
                flist = glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2]-1)+'M1[1-2].nc')+glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2])+'M0[1-9].nc')+glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2])+'M10.nc')

        # only get certain months for grid/offshore 
        elif region_name == 'grid' or region_name == 'offshore':
            if gtimep == 'fulltime':
                flist = glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2])+'M??.nc')
            if gtimep == 'springsummer':
                flist = glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2])+'M0[4-9].nc')
            if gtimep == 'junoct':
                flist = glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2])+'M0[6-9].nc')+glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2])+'M10.nc')
            if gtimep == 'sepnov':
                flist = glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2])+'M09.nc')+glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2])+'M1[0-1].nc')

            # try 2015M11-2016M10 and 2016M11-2017M10
            if gtimep == '1997':
                flist = glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2]-1)+'M1[1-2].nc')+glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2])+'M0[1-9].nc')+glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2])+'M10.nc')

        temparr = np.ones((len(flist),mask_nc.shape[0],mask_nc.shape[1]))*np.nan
        for f_i in range(len(flist)):
            temprd = np.squeeze(Dataset(flist[f_i],'r').variables[var_nc])
            cntrlm = np.squeeze(Dataset(outpath+filest+'cntrl_2012_2017_'+var_name+'_'+flist[f_i][flist[f_i].index('Y'):],'r').variables[var_nc])
            temparr[f_i] = (temprd - cntrlm)*mask_mult*s2d
            temparr[temparr==0] = np.nan

        avgp2[e_i2,y_i2] = np.nanmean(temparr) 
        stdp2[e_i2,y_i2] = np.nanstd(np.nanmean(temparr,axis=(1,2)))

    # now calculate for anth run
    if region_name == 'coast':
        if ctimep == 'fulltime':
            flist = glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2])+'M??.nc')
        if ctimep == 'springsummer': 
            flist = glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2])+'M0[4-9].nc')

        # try 2015M11-2016M10 and 2016M11-2017M10
        if ctimep == '1997': 
            flist = glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2]-1)+'M1[1-2].nc')+glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2])+'M0[1-9].nc')+glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2])+'M10.nc')

    # only get certain months for grid/offshore 
    elif region_name == 'grid' or region_name == 'offshore':
        if gtimep == 'fulltime':
            flist = glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2])+'M??.nc')
        if gtimep == 'springsummer': 
            flist = glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2])+'M0[4-9].nc')
        if gtimep == 'junoct': 
            flist = glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2])+'M0[6-9].nc')+glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2])+'M10.nc')
        if gtimep == 'sepnov': 
            flist = glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2])+'M09.nc')+glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2])+'M1[0-1].nc')

        # try 2015M11-2016M10 and 2016M11-2017M10
        if gtimep == '1997': 
            flist = glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2]-1)+'M1[1-2].nc')+glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2])+'M0[1-9].nc')+glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2])+'M10.nc')

    temparr = np.ones((len(flist),mask_nc.shape[0],mask_nc.shape[1]))*np.nan
    for f_i in range(len(flist)):
        temprd = np.squeeze(Dataset(flist[f_i],'r').variables[var_nc])
        cntrlm = np.squeeze(Dataset(outpath+filest+'cntrl_2012_2017_'+var_name+'_'+flist[f_i][flist[f_i].index('Y'):],'r').variables[var_nc])
        temparr[f_i] = (temprd - cntrlm)*mask_mult*s2d
        temparr[temparr==0] = np.nan
    avganth2[y_i2] = np.nanmean(temparr) 
    stdanth2[y_i2] = np.nanstd(np.nanmean(temparr,axis=(1,2)))

avgp98 = ((avgp1[:,0]-avganth1[0])/avganth1[0])*100
avgp99 = ((avgp1[:,1]-avganth1[1])/avganth1[1])*100
avgp16 = ((avgp2[:,0]-avganth2[0])/avganth2[0])*100
avgp17 = ((avgp2[:,1]-avganth2[1])/avganth2[1])*100

stdp98 = (stdp1[:,0]/avganth1[0])*100
stdp99 = (stdp1[:,1]/avganth1[1])*100
stdp16 = (stdp2[:,0]/avganth2[0])*100
stdp17 = (stdp2[:,1]/avganth2[1])*100

yb = -110
yt = 20

saveend = '_'+ctimep+'_'+gtimep+'.png'

# 98
figw = 16
figh = 4
axsize = 14
fig,ax = plt.subplots(1,1,figsize=[figw,figh])
ax.bar(np.arange(len(title_exp1)),avgp98,yerr=stdp98,color='white',edgecolor='k',capsize=4)
ax.set_ylim(bottom=yb,top=yt)
#if region_name == 'coast':
#    ax.set_ylim(bottom=-0.5,top=14)
#if region_name == 'grid':
#    ax.set_ylim(bottom=-1.5,top=4.7)
#if region_name == 'offshore':
#    ax.set_ylim(bottom=-1.5,top=4.7)
ax.plot(range(len(title_exp1)),np.ones((len(title_exp1)))*0)
ax.set_xticks(range(len(title_exp1)))
ax.set_xticklabels(title_exp1,fontsize=axsize)
ax.set_ylabel('% change integrated NPP',fontsize=axsize)
ax.set_title('1998',fontsize=axsize)
ax.tick_params(axis='both',which='major',labelsize=axsize)
savename = var_name+'_avg_4years_byyear_perc'+region_name+'_98'+saveend+'.png'
fig.savefig(savepath+savename,bbox_inches='tight')
print(savename)

# 99
figw = 16
figh = 4
axsize = 14
fig,ax = plt.subplots(1,1,figsize=[figw,figh])
ax.bar(np.arange(len(title_exp1)),avgp99,yerr=stdp99,color='white',edgecolor='k',capsize=4)
ax.set_ylim(bottom=yb,top=yt)
#if region_name == 'coast':
#    ax.set_ylim(bottom=-0.5,top=14)
#if region_name == 'grid':
#    ax.set_ylim(bottom=-1.5,top=4.7)
#if region_name == 'offshore':
#    ax.set_ylim(bottom=-1.5,top=4.7)
ax.plot(range(len(title_exp1)),np.ones((len(title_exp1)))*0)
ax.set_xticks(range(len(title_exp1)))
ax.set_xticklabels(title_exp1,fontsize=axsize)
ax.set_ylabel('% change integrated NPP',fontsize=axsize)
ax.set_title('1999',fontsize=axsize)
ax.tick_params(axis='both',which='major',labelsize=axsize)
savename = var_name+'_avg_4years_byyear_perc'+region_name+'_99'+saveend+'.png'
fig.savefig(savepath+savename,bbox_inches='tight')
print(savename)

# 16
figw = 12
figh = 4
axsize = 14
fig,ax = plt.subplots(1,1,figsize=[figw,figh])
ax.bar(range(len(title_exp2)),avgp16,yerr=stdp16,color='gray',capsize=4)
ax.set_ylim(bottom=yb,top=yt)
#if region_name == 'coast':
#    ax.set_ylim(bottom=-0.5,top=14)
#if region_name == 'grid':
#    ax.set_ylim(bottom=-1.5,top=4.7)
#if region_name == 'offshore':
#    ax.set_ylim(bottom=-1.5,top=4.7)
ax.plot(range(len(title_exp2)),np.ones((len(title_exp2)))*0)
ax.set_xticks(range(len(title_exp2)))
ax.set_xticklabels(title_exp2,fontsize=axsize)
ax.set_ylabel('% change integrated NPP',fontsize=axsize)
ax.set_title('2016',fontsize=axsize)
ax.tick_params(axis='both',which='major',labelsize=axsize)
savename = var_name+'_avg_4years_byyear_perc'+region_name+'_16'+saveend+'.png'
fig.savefig(savepath+savename,bbox_inches='tight')
print(savename)

# 17
figw = 12
figh = 4
axsize = 14
fig,ax = plt.subplots(1,1,figsize=[figw,figh])
ax.bar(range(len(title_exp2)),avgp17,yerr=stdp17,color='gray',capsize=4)
ax.set_ylim(bottom=yb,top=yt)
#if region_name == 'coast':
#    ax.set_ylim(bottom=-0.5,top=14)
#if region_name == 'grid':
#    ax.set_ylim(bottom=-1.5,top=4.7)
#if region_name == 'offshore':
#    ax.set_ylim(bottom=-1.5,top=4.7)
ax.plot(range(len(title_exp2)),np.ones((len(title_exp2)))*0)
ax.set_xticks(range(len(title_exp2)))
ax.set_xticklabels(title_exp2,fontsize=axsize)
ax.set_ylabel('% change integrated NPP',fontsize=axsize)
ax.set_title('2017',fontsize=axsize)
ax.tick_params(axis='both',which='major',labelsize=axsize)
savename = var_name+'_avg_4years_byyear_perc'+region_name+'_17'+saveend+'.png'
fig.savefig(savepath+savename,bbox_inches='tight')
print(savename)
