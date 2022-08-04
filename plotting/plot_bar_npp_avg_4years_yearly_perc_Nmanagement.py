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

# months to average over for offshore/full bight
# 7 - 11 (hard coded in glob loop)

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
       #'pndn50', 
       #'pndn90', 
       'FNDN_only']
       #'fndn50',
       #'fndn90']

exp2 = ['PNDN_only_realistic',
       #'pndn50_realistic', 
       #'pndn90_realistic',
       'FNDN_only_realistic']

title_exp = ['50% N Red.',
             '',
             #'50% N Red.\n50% Recy.',
             #'',
             #'50% N Red.\n90% Recy.',
             #'',
             '85% N Red.',
             '']
             #'85% N Red.\n50% Recy.',
             #'85% N Red.\n90% Recy.']

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
    regtitle = 'Offshore'

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

fileen = '_'+var_name+'.nc'

# 1998-1999
avganth1 = np.ones((len(yearlist1)))*np.nan
avgp1 = np.ones((len(exp1),len(yearlist1)))*np.nan
for y_i1 in range(len(yearlist1)):
    for e_i1 in range(len(exp1)):
        # get all months
        if region_name == 'coast':
            flist = glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1])+'M??.nc')
        # only get certain months for grid/offshore 
        elif region_name == 'grid' or region_name == 'offshore':
            flist = glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1])+'M??.nc')
            #flist = glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1])+'M0[6-9].nc')+glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1])+'M10.nc')
            #flist = glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1])+'M09.nc')+glob.glob(outpath+filest+exp1[e_i1]+'_'+var_name+'_Y'+str(yearlist1[y_i1])+'M1[0-1].nc')
        temparr = np.ones((len(flist),mask_nc.shape[0],mask_nc.shape[1]))*np.nan
        for f_i in range(len(flist)):
            temprd = np.squeeze(Dataset(flist[f_i],'r').variables[var_nc])
            cntrlm = np.squeeze(Dataset(outpath+filest+'cntrl_initap_'+var_name+'_'+flist[f_i][flist[f_i].index('Y'):],'r').variables[var_nc])
            temparr[f_i] = (temprd - cntrlm)*mask_mult*s2d
            temparr[temparr==0] = np.nan

        avgp1[e_i1,y_i1] = np.nanmean(temparr) 

    # now calculate for anth run
    if region_name == 'coast':
        flist = glob.glob(outpath+filest+'loads1617_'+var_name+'_Y'+str(yearlist1[y_i1])+'M??.nc')
    # only get certain months for grid/offshore 
    elif region_name == 'grid' or region_name == 'offshore':
        flist = glob.glob(outpath+filest+'loads1617_'+var_name+'_Y'+str(yearlist1[y_i1])+'M??.nc')
        #flist = glob.glob(outpath+filest+'loads1617_'+var_name+'_Y'+str(yearlist1[y_i1])+'M0[6-9].nc')+glob.glob(outpath+filest+'loads1617_'+var_name+'_Y'+str(yearlist1[y_i1])+'M10.nc')
        #flist = glob.glob(outpath+filest+'loads1617_'+var_name+'_Y'+str(yearlist1[y_i1])+'M09.nc')+glob.glob(outpath+filest+'loads1617_'+var_name+'_Y'+str(yearlist1[y_i1])+'M1[0-1].nc')
    temparr = np.ones((len(flist),mask_nc.shape[0],mask_nc.shape[1]))*np.nan
    for f_i in range(len(flist)):
        temprd = np.squeeze(Dataset(flist[f_i],'r').variables[var_nc])
        cntrlm = np.squeeze(Dataset(outpath+filest+'cntrl_initap_'+var_name+'_'+flist[f_i][flist[f_i].index('Y'):],'r').variables[var_nc])
        temparr[f_i] = (temprd - cntrlm)*mask_mult*s2d
        temparr[temparr==0] = np.nan
    avganth1[y_i1] = np.nanmean(temparr) 
    

# 2016-2017
avganth2 = np.ones((len(yearlist2)))*np.nan
avgp2 = np.ones((len(exp2),len(yearlist2)))*np.nan
for y_i2 in range(len(yearlist2)):
    for e_i2 in range(len(exp2)):
        if region_name == 'coast':
            flist = glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2])+'M??.nc')
        # only get certain months for grid/offshore 
        elif region_name == 'grid' or region_name == 'offshore':
            flist = glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2])+'M??.nc')
            #flist = glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2])+'M0[6-9].nc')+glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2])+'M10.nc')
            #flist = glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2])+'M09.nc')+glob.glob(outpath+filest+exp2[e_i2]+'_'+var_name+'_Y'+str(yearlist2[y_i2])+'M1[0-1].nc')
        temparr = np.ones((len(flist),mask_nc.shape[0],mask_nc.shape[1]))*np.nan
        for f_i in range(len(flist)):
            temprd = np.squeeze(Dataset(flist[f_i],'r').variables[var_nc])
            cntrlm = np.squeeze(Dataset(outpath+filest+'cntrl_2012_2017_'+var_name+'_'+flist[f_i][flist[f_i].index('Y'):],'r').variables[var_nc])
            temparr[f_i] = (temprd - cntrlm)*mask_mult*s2d
            temparr[temparr==0] = np.nan

        avgp2[e_i2,y_i2] = np.nanmean(temparr) 

    # now calculate for anth run
    if region_name == 'coast':
        flist = glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2])+'M??.nc')
    # only get certain months for grid/offshore 
    elif region_name == 'grid' or region_name == 'offshore':
        flist = glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2])+'M??.nc')
        #flist = glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2])+'M0[6-9].nc')+glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2])+'M10.nc')
        #flist = glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2])+'M09.nc')+glob.glob(outpath+filest+'fulll_2012_2017_'+var_name+'_Y'+str(yearlist2[y_i2])+'M1[0-1].nc')
    temparr = np.ones((len(flist),mask_nc.shape[0],mask_nc.shape[1]))*np.nan
    for f_i in range(len(flist)):
        temprd = np.squeeze(Dataset(flist[f_i],'r').variables[var_nc])
        cntrlm = np.squeeze(Dataset(outpath+filest+'cntrl_2012_2017_'+var_name+'_'+flist[f_i][flist[f_i].index('Y'):],'r').variables[var_nc])
        temparr[f_i] = (temprd - cntrlm)*mask_mult*s2d
        temparr[temparr==0] = np.nan
    avganth2[y_i2] = np.nanmean(temparr) 

avgp = np.ones((len(title_exp)))*np.nan
stdp = np.ones((len(title_exp)))*np.nan

avganth = np.ones((len(title_exp)))*np.nan
stdanth = np.ones((len(title_exp)))*np.nan

avgp[0] = np.nanmean((avgp1[0]))
avgp[1] = np.nanmean((avgp2[0]))
avgp[2] = np.nanmean((avgp1[1]))
avgp[3] = np.nanmean((avgp2[1]))
#avgp[4] = np.nanmean((avgp1[2]))
#avgp[5] = np.nanmean((avgp2[2]))
#avgp[6] = np.nanmean((avgp1[3]))
#avgp[7] = np.nanmean((avgp2[3]))
#avgp[8] = np.nanmean((avgp1[4]))
#avgp[9] = np.nanmean((avgp1[5]))

stdp[0] = np.nanstd((avgp1[0]))
stdp[1] = np.nanstd((avgp2[0]))
stdp[2] = np.nanstd((avgp1[1]))
stdp[3] = np.nanstd((avgp2[1]))
#stdp[4] = np.nanstd((avgp1[2]))
#stdp[5] = np.nanstd((avgp2[2]))
#stdp[6] = np.nanstd((avgp1[3]))
#stdp[7] = np.nanstd((avgp2[3]))
#stdp[8] = np.nanstd((avgp1[4]))
#stdp[9] = np.nanstd((avgp1[5]))

ind = [1,3]

avganth[:] = np.nanmean((avganth1))
avganth[ind] = np.nanmean((avganth2))

stdanth[:] = np.nanstd((avganth1))
stdanth[ind] = np.nanstd((avganth2))

avgperc = ((avgp-avganth)/avganth)*100
#stdperc = ((stdp-avganth)/avganth)*100
#stdperc = (stdp/avgp)*100
stdperc = (stdp/avganth)*100

figw = 10
figh = 4

axsize = 14

#indices = [0,2,4,6,8,9]
indices = [0,2]

fig,ax = plt.subplots(1,1,figsize=[figw,figh])
ax.bar(np.arange(len(title_exp))[indices],avgperc[indices],yerr=stdperc[indices],color='white',edgecolor='k',capsize=4)
ax.bar(np.arange(len(title_exp))[ind],avgperc[ind],yerr=stdperc[ind],color='gray',capsize=4)

ax.set_ylim(bottom=-120,top=10)
#if region_name == 'coast':
#    ax.set_ylim(bottom=-0.5,top=14)
#if region_name == 'grid':
#    ax.set_ylim(bottom=-1.5,top=4.7)
#if region_name == 'offshore':
#    ax.set_ylim(bottom=-1.5,top=4.7)

ax.plot(range(len(title_exp)),np.ones((len(title_exp)))*0)

xloc = np.arange(len(title_exp))+.5
#xloc[-2:] = xloc[-2:]-0.5

ax.set_xticks(xloc)
ax.set_xticklabels(title_exp,fontsize=axsize)
ax.set_ylabel('Integrated NPP 100 m '+cblabel,fontsize=axsize)
ax.set_title(regtitle,fontsize=axsize)
ax.tick_params(axis='both',which='major',labelsize=axsize)

savename = var_name+'_avg_4years_yearly_perc'+region_name+'_Nmanagement.png'
fig.savefig(savepath+savename,bbox_inches='tight')
print(savename)
#plt.close()

