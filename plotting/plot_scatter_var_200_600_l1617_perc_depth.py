# sum up all negative values in a mask
# used to see impact of treatment level on oxygen
# make my own 5-10 km mask circles?
import sys
import os
sys.path.append('/data/project3/minnaho/global/')
import l2grid as l2grid
import numpy as np
from netCDF4 import Dataset,num2date
import matplotlib.pyplot as plt

# plot fresh vs nutrients vs control vs full
plt.ion()

savepath = './figs/scatter/'
region_name = 'grid'
# index of depths: 0 = 2 m, 1 = 4 m, etc
depl0 = 50  # choose depth layer to start from, actual depth *2
depl1 = 110 # choose depth layer to end, actual depth *2
            # but if >100, then ((x-100)*40)+200 is the depth

# ROMS output location
#outpath = '/data/project6/minnaho/opc_scenarios/ext_depth_200/'
outpath = '/data/project6/minnaho/opc_scenarios/ext_depth_200_monthly/'
outpath1 = '/data/project6/minnaho/opc_scenarios/ext_depth_600/'

# roms var
var_name = 'O2' 
var_nc = 'var' 

#year_month = 'spring1998'
year_month = 'fullts'

#filest = 'ext_0_200_'+var_name+'_'+year_month+'_'
filest = 'avg_fullts_0_200_'+var_name+'_'
fileen = '.nc'
filest1 = 'avg_fullts_200_600_'+var_name+'_'
fileen1 = '.nc'

# scenario names 
#exp = ['PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
#title_exp = ['PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']
#title_exp = ['50% N\nReduction','50% N\nReduction\n50% Recycle','50% N\nReduction\n90% Recycle','85% N\nReduction','85% N\nReduction\n50% Recycle','85% N\nReduction\n90% Recycle',]
exp = ['PNDN_only','FNDN_only']
title_exp = ['50% N\nReduction','85% N\nReduction']

mmolm3_to_mgl = 32./1000

# region masks
mask_nc = l2grid.mask_nc
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
    mask_mult = mask_nc
    regtitle = 'full SCB'

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

pm_nc = l2grid.pm_nc
pn_nc = l2grid.pn_nc

mask_nc[mask_nc==0] = np.nan # get rid of land volumes
xisize = (1/pm_nc)*mask_nc
etasize = (1/pn_nc)*mask_nc

fpath = []
for e_i in range(len(exp)):
    fpath.append(outpath+filest+exp[e_i]+fileen)

fpath1 = []
for e_i in range(len(exp)):
    fpath1.append(outpath1+filest1+exp[e_i]+fileen)

# read in control to subtract from
if depl1 <= 100:
    roms_cnt = (np.array(Dataset(outpath+filest+'cntrl'+fileen,'r').variables[var_nc])*mask_mult)[:,depl0:depl1,:,:]
    roms_cnt_vert = np.nanmean(roms_cnt,axis=1)
    roms_cnt_mean = np.nanmean(roms_cnt_vert)
    cntrl_mean = roms_cnt_mean*mmolm3_to_mgl

else:
    roms_cnt = (np.array(Dataset(outpath+filest+'cntrl'+fileen,'r').variables[var_nc])*mask_mult)[:,depl0:,:,:]
    roms_cnt1 = (np.array(Dataset(outpath1+filest1+'cntrl'+fileen,'r').variables[var_nc])*mask_mult)[:,:(depl1-100),:,:]
    #roms_cnt1_vert = np.nanmean(roms_cnt1,axis=1)
    #roms_cnt1_mean = np.nanmean(roms_cnt1_vert)
    
    #cntrl_mean = np.nanmean((roms_cnt_mean,roms_cnt1_mean))*mmolm3_to_mgl
    
    # concat two files
    cnt_concat = np.concatenate((roms_cnt,roms_cnt1),axis=1)
    cnt_norm = np.nansum(cnt_concat,axis=1)/cnt_concat.shape[1]
    cntrl_mean = np.nanmean(cnt_norm)*mmolm3_to_mgl

# l1617 upper limit
if depl1 <= 100:
    roms_l16 = (np.array(Dataset(outpath+filest+'l1617'+fileen,'r').variables[var_nc])*mask_mult)[:,depl0:depl1,:,:]
    roms_l16_vert = np.nanmean(roms_l16,axis=1)
    roms_l16_mean = np.nanmean(roms_l16_vert)
    l1617_mean = roms_l16_mean*mmolm3_to_mgl
    
else:
    roms_l16 = (np.array(Dataset(outpath+filest+'l1617'+fileen,'r').variables[var_nc])*mask_mult)[:,depl0:,:,:]
    roms_l161 = (np.array(Dataset(outpath1+filest1+'l1617'+fileen,'r').variables[var_nc])*mask_mult)[:,:depl1-100,:,:]
    roms_l161_vert = np.nanmean(roms_l161,axis=1)
    roms_l161_mean = np.nanmean(roms_l161_vert)

    #l1617_mean = np.nanmean((roms_l16_mean,roms_l161_mean))*mmolm3_to_mgl
    
    # concat two files
    l16_concat = np.concatenate((roms_l16,roms_l161),axis=1)
    l16_norm = np.nansum(l16_concat,axis=1)/l16_concat.shape[1]
    l1617_mean = np.nanmean(l16_norm)*mmolm3_to_mgl
    
l1617_perc = ((l1617_mean-cntrl_mean)/cntrl_mean)*100

#figw = 14
figw = 10
figh = 4

axis_tick_size = 14

#savename = 'sum_posneg_'+var_name+'_'+year_month+'_'+region_name+'_0_600_'+exp[-1]+'_bar_norm_'+str(depl0*2)+'_600m.png'
if depl1 <= 100:
    savename = 'sum_posneg_'+var_name+'_'+year_month+'_'+region_name+'_0_600_'+exp[-1]+'_bar_norm_'+str(depl0*2)+'_'+str(depl1*2)+'.png'
else:
    savename = 'sum_posneg_'+var_name+'_'+year_month+'_'+region_name+'_0_600_'+exp[-1]+'_bar_norm_'+str(depl0*2)+'_'+str(((depl1-100)*40)+200)+'.png'
#savename = 'sum_neg_'+var_name+'_'+year_month+'_'+region_name+'_200_'+exp[-1]+'.png'
fig,ax = plt.subplots(1,1,figsize=[figw,figh])


for n_i in range(len(exp)):

    #0-200

    if depl1 <= 100:
        roms_var_read = (np.array(Dataset(fpath[n_i],'r').variables[var_nc])*mask_mult)[:,depl0:depl1,:,:]
    #200-600
    #roms_var_read1 = np.array(Dataset(fpath1[n_i],'r').variables[var_nc])*mask_mult

        roms_neg_vert = np.nanmean(roms_var_read,axis=1)
        roms_neg = np.nanmean(roms_neg_vert)*mmolm3_to_mgl

    #roms_var = roms_var_read - roms_cnt
    #roms_var1 = roms_var_read1 - roms_cnt1

    else:
        # 0-200
        roms_var_read = (np.array(Dataset(fpath[n_i],'r').variables[var_nc])*mask_mult)[:,depl0:,:,:]
        #200-600
        roms_var_read1 = (np.array(Dataset(fpath1[n_i],'r').variables[var_nc])*mask_mult)[:,:depl1-100,:,:]
        # concat two files
        roms_concat = np.concatenate((roms_var_read,roms_var_read1),axis=1)
        roms_var_norm = np.nansum(roms_concat,axis=1)/roms_concat.shape[1]
        roms_neg = np.nanmean(roms_var_norm)*mmolm3_to_mgl
        print('maximum increase in oxygen ',exp[e_i],str(np.nanmax(((roms_concat-l16_concat)/l16_concat)*100)+'%'))

    # take average
    #roms_neg_vert = np.nanmean(roms_var_read,axis=1)
    #roms_neg_mean = np.nanmean(roms_neg_vert)

    #roms_neg_vert1 = np.nanmean(roms_var_read1,axis=1)
    #roms_neg_mean1 = np.nanmean(roms_neg_vert1)

    #roms_neg = np.nanmean((roms_neg_mean,roms_neg_mean1))*mmolm3_to_mgl

    # normalize by number of depth slices
    #roms_neg0 = np.nansum(roms_var)/roms_var.shape[1]
    #roms_neg1 = np.nansum(roms_var1)/roms_var1.shape[1]

    #roms_neg_sum = roms_neg0+roms_neg1

    # divide by bight area to get volume
    #roms_neg_sum = (roms_neg0+roms_neg1)/np.nansum(xisize*etasize)

    #if exp[n_i] == 'l1617':
    #    roms_l = roms_neg_sum
    #    print(exp[n_i]+': ',str(roms_l))
    #else: 
    #    print(exp[n_i]+': ',str(roms_neg_sum))
    #    #roms_neg = ((roms_neg_sum-roms_l)/roms_l)*100
    #    roms_neg = roms_neg_sum-roms_l
    #    ax.scatter(n_i,roms_neg)

    print(exp[n_i],str(roms_neg))
    #ax.bar(n_i,roms_neg-l1617_mean,bottom=l1617_mean)
    ax.bar(n_i,((roms_neg-cntrl_mean)/cntrl_mean)*100)

#ax.set_yscale('log')
#ax.set_ybound(lower=3E5,upper=5E7)
#ax.set_ylim(bottom=l1617_mean-(l1617_mean*0.01),top=cntrl_mean+(cntrl_mean*0.01))
#ax.set_ylim(bottom=l1617_perc-(l1617_perc*0.01),top=-3)
ax.set_ylim(bottom=-5,top=-4.2)
ax.set_xticks(range(len(exp)))
#ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*cntrl_mean,linestyle='--',color='purple')
ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*l1617_perc,linestyle='--',color='green')

ax.set_xticklabels(title_exp)
ax.set_xlabel('Scenario',fontsize=axis_tick_size)
#ax.set_ylabel('Sum of O2 change mmol m$^{-2}$',fontsize=axis_tick_size)
#ax.set_ylabel('Average O2 change mg L$^{-1}$',fontsize=axis_tick_size)
ax.set_ylabel('% O2 change',fontsize=axis_tick_size)
ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)
ax.tick_params(axis='both',which='minor',labelsize=axis_tick_size)

#fig.suptitle('Sum of Negative O2 '+year_month+' Average '+regtitle,fontsize=axis_tick_size)


fig.savefig(savepath+savename,bbox_inches='tight')
print(savename)
#plt.close()

