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

savepath = './figs/2years/'
region_name = 'grid'

# coast time period
ctimep = '1997'
# grid/offshore time period
gtimep = 'fulltime'

# ROMS output location
outpath = '/data/project6/minnaho/opc_scenarios/bgc_flux/'

# roms var
var_name = 'npp' 
var_nc = 'var_int' 
cblabel = 'mmol m$^{-2}$ d$^{-1}$'

filest = 'int_avg_100m_50m_'

yearlist1 = [1998,1999]
yearlist2 = [2016,2017]

exp1 = ['cntrl_initap_realistic','fulll_2012_2017']

exp2 = ['PNDN_only_realistic',
       'pndn50_fixriver', 
       'pndn90_fixriver',
       'FNDN_only_realistic',
        'fndn50_fixriver',
        'fndn90_fixriver']

title_exp = ['50% N Red.',
             '50% N Red.\n50% Recy.',
             '50% N Red.\n90% Recy.',
             '85% N Red.',
             '85% N Red.\n50% Recy.',
             '85% N Red.\n90% Recy.'
             ]

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
    '''
    mask7[mask9==1] = 1
    mask_mult = mask7
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

# 2016-2017
avganth2 = np.ones((len(yearlist2)))*np.nan
avgp2 = np.ones((len(exp2),len(yearlist2)))*np.nan
stdanth2 = np.ones((len(yearlist2)))*np.nan
stdp2 = np.ones((len(exp2)+2,len(yearlist2)))*np.nan

# raw npp values for cntrl/anth
npp_maps_cntrl = np.ones((24,mask_nc.shape[0],mask_nc.shape[1])) 
npp_maps_anth = np.ones((24,mask_nc.shape[0],mask_nc.shape[1])) 

# raw npp values for each scenario
npp_maps = np.ones((len(exp2),24,mask_nc.shape[0],mask_nc.shape[1]))*np.nan

# scenario-anth npp
npp_diff_anth = np.ones((len(exp2),24,mask_nc.shape[0],mask_nc.shape[1]))*np.nan

# scenario-cntrl npp
npp_diff_cntrl = np.ones((len(exp2),24,mask_nc.shape[0],mask_nc.shape[1]))*np.nan

# monthly means: raw, scenario-anth, scenario-cntrl
npp_monthly = np.ones((len(exp2),24))*np.nan
npp_monthly_anth = np.ones((len(exp2),24))*np.nan
npp_monthly_cntrl = np.ones((len(exp2),24))*np.nan

for e_i in range(len(exp1)):
    flist = glob.glob(outpath+filest+exp1[e_i]+'_'+var_name+'_Y*')
    for f_i in range(len(flist)):
        if e_i == 0:
            npp_maps_cntrl[f_i] = np.squeeze(Dataset(flist[f_i],'r')[var_nc])
        if e_i == 1:
            npp_maps_anth[f_i] = np.squeeze(Dataset(flist[f_i],'r')[var_nc])

npp_mean_monthly_cntrl = np.nanmean(npp_maps_cntrl*mask_mult,axis=(1,2))*s2d

npp_mean_anth_cntrl = np.nanmean((npp_maps_anth-npp_maps_cntrl)*mask_mult)*s2d

npp_mean_anth = np.nanmean(npp_maps_anth*mask_mult)*s2d
npp_mean_monthly_anth = np.nanmean(npp_maps_anth*mask_mult,axis=(1,2))*s2d

for e_i in range(len(exp2)):
    flist = glob.glob(outpath+filest+exp2[e_i]+'_'+var_name+'_Y*')
    for f_i in range(len(flist)):
        npp_maps[e_i,f_i] = np.squeeze(Dataset(flist[f_i],'r')[var_nc])
        npp_diff_anth[e_i,f_i] = npp_maps[e_i,f_i]-npp_maps_anth[f_i]
        npp_diff_cntrl[e_i,f_i] = npp_maps[e_i,f_i]-npp_maps_cntrl[f_i]
    # monthly mean npp
    npp_monthly[e_i] = np.nanmean(npp_maps[e_i,:,:,:]*mask_mult,axis=(1,2))*s2d
    npp_monthly_anth[e_i] = np.nanmean(npp_diff_anth[e_i,:,:,:]*mask_mult,axis=(1,2))
    npp_monthly_cntrl[e_i] = np.nanmean(npp_diff_cntrl[e_i,:,:,:]*mask_mult,axis=(1,2))

npp_mean_sce_anth = np.nanmean(npp_diff_anth,axis=(1,2,3))*s2d

npp_mean = np.nanmean(npp_monthly,axis=1)
npp_monthly_05 = np.nanpercentile(npp_monthly,5,axis=1)
npp_monthly_95 = np.nanpercentile(npp_monthly,95,axis=1)

npp_05_sce_anth = np.nanpercentile(npp_diff_anth,5,axis=(1,2,3))
npp_95_sce_anth = np.nanpercentile(npp_diff_anth,95,axis=(1,2,3))

npp_2016 = np.nanmean(npp_monthly[:,:12],axis=1)
npp_2017 = np.nanmean(npp_monthly[:,12:],axis=1)

npp25 = np.nanpercentile(np.array((npp_2016,npp_2017)),25,axis=0)
npp75 = np.nanpercentile(np.array((npp_2016,npp_2017)),75,axis=0)

# subtract (scenario-anth)-(anth-cntrl)/(anth-cntrl)
perc = ((npp_mean-npp_mean_anth)/npp_mean_anth_cntrl)*100 
perc25 = ((npp25-npp_mean_anth)/npp_mean_anth_cntrl)*100 
perc75 = ((npp75-npp_mean_anth)/npp_mean_anth_cntrl)*100 

plt.ion()

figw = 12
figh = 4
axfont = 14

nman_name = ['0% N Red.','50% N Red.','90% N Red.']

fig,ax = plt.subplots(1,1,figsize=[figw,figh])
ax.errorbar([0,50,90],[0,perc[0],perc[3]],yerr=np.array((np.abs((0,perc25[0],perc25[3])),np.abs((0,perc75[0],perc75[3])))))
ax.scatter([0,50,90],[0,perc[0],perc[3]])
ax.plot([0,100],np.repeat(-100,2),color='k',linestyle='--')
ax.plot([0,100],np.repeat(0,2),color='k',linestyle='--')
ax.set_ylim([-105,5])
ax.set_ylabel('% change in NPP from ANTH',fontsize=axfont)
ax.set_xlabel('% Inorganic N Reduction',fontsize=axfont)
#ax.set_xticklabels(nman_name,fontsize=axfont)
ax.tick_params(axis='both',which='major',right=True,labelsize=axfont)

dates = ['11/2015','12/2015','1/2016','2/2016','3/2016','4/2016','5/2016','6/2016','7/2016','8/2016','9/2016','10/2016','11/2016','12/2016','1/2017','2/2017','3/2017','4/2017','5/2017','6/2017','7/2017','8/2017','9/2017','10/2017']

fig,ax = plt.subplots(6,1,figsize=[7,16],sharex=True)
for e_i in range(len(npp_monthly)):
    #ax.flat[e_i].plot(dates,npp_monthly[e_i])
    ax.flat[e_i].plot(dates,npp_monthly_anth[e_i]*s2d)

ax.flat[5].set_xticklabels(dates,rotation=60,ha='right')
for a_i in range(len(ax.flat)):
    ax.flat[e_i].set_ylim([-65,50])


'''
#figw = 16
figw = 12
figh = 4

axsize = 14

ind1 = list(range(0,6))
ind2 = list(range(6,len(title_exp)))

nman_name = ['No N Red.','50% N Red.','85% N Red.']

fig,ax = plt.subplots(1,1,figsize=[figw,figh])
#ax.plot(nman_name,pnmanavg,color='k',linewidth=3)
ax.errorbar(nman_name,pnmanavg,yerr=pnmanstd,linewidth=2)
#ax.bar(np.arange(len(title_exp))[ind2],avgnum2,yerr=stdnum2,color='gray',capsize=4)

ax.plot(np.arange(len(nman_name)),np.zeros((len(nman_name))),color='k',linestyle='--')
ax.plot(np.arange(len(nman_name)),np.ones((len(nman_name)))*-100,color='k',linestyle='--')

#if region_name == 'coast':
#    ax.set_ylim(bottom=-0.5,top=14)
#if region_name == 'grid':
#    ax.set_ylim(bottom=-1.5,top=3.5)
#if region_name == 'offshore':
#    ax.set_ylim(bottom=-1.5,top=3.5)

ax.set_ylim(bottom=-115,top=10)
#ax.yaxis.set_ticks(np.arange(-100,20,20))

#ax.plot(np.arange(len(title_exp))[ind1],np.ones((len(title_exp)))[ind1]*avganthnum1)
#ax.fill_between(np.arange(len(title_exp))[ind1],np.ones((len(title_exp)))[ind1]*(avganthnum1+stdanthnum1),np.ones((len(title_exp)))[ind1]*(avganthnum1-stdanthnum1),alpha=0.3)
#
#ax.plot(np.arange(len(title_exp))[ind2],np.ones((len(title_exp)))[ind2]*avganthnum2)
#ax.fill_between(np.arange(len(title_exp))[ind2],np.ones((len(title_exp)))[ind2]*(avganthnum2+stdanthnum2),np.ones((len(title_exp)))[ind2]*(avganthnum2-stdanthnum2),alpha=0.3)

ax.set_xticklabels(nman_name,fontsize=axsize)
ax.set_ylabel('% Change Algal Production',fontsize=axsize)
ax.set_title(regtitle,fontsize=axsize)
ax.tick_params(axis='both',which='major',right=True,labelsize=axsize)

savename = 'line_'+var_name+'_avg_n4_perc_'+region_name+'_'+ctimep+'_'+gtimep+'_rearrange.png'

fig.savefig(savepath+savename,bbox_inches='tight')
print(savename)
#plt.close()
'''

