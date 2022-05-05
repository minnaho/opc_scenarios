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

perc = False

region_file = ['grid','mask1','mask2','mask3','mask4','mask5','mask6','mask7','mask8','mask9']
regtitle = ['SCB','South SD','Nth SD & Sth OC','Los Angeles coast','North coast','Channel Island','S. Barbara B','S. Monica/Pedro B','Islands Basins','SD Terrasse']

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

figw = 8
figh = 12

axis_tick_size = 14

if depl1 <= 100:
    if perc == True:
        savename = 'barh_regions_change_'+var_name+'_'+year_month+'_0_600_'+exp[-1]+'_'+str(depl0*2)+'_'+str(depl1*2)+'_perc.png'
    else:
        savename = 'barh_regions_change_'+var_name+'_'+year_month+'_0_600_'+exp[-1]+'_'+str(depl0*2)+'_'+str(depl1*2)+'.png'
else:
    if perc == True:
        savename = 'barh_regions_change_'+var_name+'_'+year_month+'_0_600_'+exp[-1]+'_'+str(depl0*2)+'_'+str(((depl1-100)*40)+200)+'_perc.png'
    else:
        savename = 'barh_regions_change_'+var_name+'_'+year_month+'_0_600_'+exp[-1]+'_'+str(depl0*2)+'_'+str(((depl1-100)*40)+200)+'.png'

fig,ax = plt.subplots(1,1,figsize=[figw,figh])



pm_nc = l2grid.pm_nc
pn_nc = l2grid.pn_nc

mask_nc[mask_nc==0] = np.nan # get rid of land volumes
xisize = (1/pm_nc)*mask_nc
etasize = (1/pn_nc)*mask_nc

pltval = np.ones((len(region_file),len(exp)+1))*np.nan

fpath = []
for e_i in range(len(exp)):
    fpath.append(outpath+filest+exp[e_i]+fileen)

fpath1 = []
for e_i in range(len(exp)):
    fpath1.append(outpath1+filest1+exp[e_i]+fileen)

for r_i in range(len(region_file)):
    region_name = region_file[r_i]
    if region_name == 'grid':
        mask_mult = mask_nc
        regtitle = 'full SCB'
    if region_name == 'mask0':
        mask_mult = mask0
        regtitle = 'Tijuana R'
    if region_name == 'mask1':
        mask_mult = mask1
        regtitle = 'South SD'
    if region_name == 'mask2':
        mask_mult = mask2
        regtitle = 'Nth SD & Sth OC'
    if region_name == 'mask3':
        mask_mult = mask3
        regtitle = 'Los Angeles coast'
    if region_name == 'mask4':
        mask_mult = mask4
        regtitle = 'North coast'
    if region_name == 'mask5':
        mask_mult = mask5
        regtitle = 'Channel Island'
    if region_name == 'mask6':
        mask_mult = mask6
        regtitle = 'S. Barbara B'
    if region_name == 'mask7':
        mask_mult = mask7
        regtitle = 'S. Monica/Pedro B'
    if region_name == 'mask8':
        mask_mult = mask8
        regtitle = 'Islands Basins'
    if region_name == 'mask9':
        mask_mult = mask9
        regtitle = 'SD Terrasse'

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
    if perc == True:
        pltval[r_i,0] = l1617_perc
    else:
        pltval[r_i,0] = l1617_mean - cntrl_mean

    for n_i in range(len(exp)):
    
        #0-200
        if depl1 <= 100:
            roms_var_read = (np.array(Dataset(fpath[n_i],'r').variables[var_nc])*mask_mult)[:,depl0:depl1,:,:]
            roms_neg_vert = np.nanmean(roms_var_read,axis=1)
            roms_neg = np.nanmean(roms_neg_vert)*mmolm3_to_mgl
    
        else:
            # 0-200
            roms_var_read = (np.array(Dataset(fpath[n_i],'r').variables[var_nc])*mask_mult)[:,depl0:,:,:]
            #200-600
            roms_var_read1 = (np.array(Dataset(fpath1[n_i],'r').variables[var_nc])*mask_mult)[:,:depl1-100,:,:]
            # concat two files
            roms_concat = np.concatenate((roms_var_read,roms_var_read1),axis=1)
            roms_var_norm = np.nansum(roms_concat,axis=1)/roms_concat.shape[1]
            roms_neg = np.nanmean(roms_var_norm)*mmolm3_to_mgl

        # +1 because no reduction/loads 1617 is first index
        if perc == True:
            pltval[r_i,n_i+1] = ((roms_neg - cntrl_mean)/cntrl_mean)*100
        else:
            pltval[r_i,n_i+1] = roms_neg - cntrl_mean

        print(exp[n_i],str(roms_neg))

width = 0.1
y_pos = np.arange(len(region_file))
ax.barh(y_pos,pltval[:,0],color='yellow',label='No reduction')
ax.barh(y_pos+width,pltval[:,1],color='red',label='50% reduction')
ax.barh(y_pos+(2*width),pltval[:,2],color='blue',label='85% reduction')

#ax.set_ylim(bottom=-5,top=-4.2)
ax.set_yticks(y_pos+width/3)

ax.set_yticklabels(regtitle)
if perc == True:
    ax.set_xlabel('% O2 change',fontsize=axis_tick_size)
else:
    ax.set_xlabel('O2 change mg L$^{-1}$',fontsize=axis_tick_size)

ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)
ax.tick_params(axis='both',which='minor',labelsize=axis_tick_size)

fig.savefig(savepath+savename,bbox_inches='tight')
print(savename)
#plt.close()

