# bar plot of percent diff from fulll-cntrl
# -100% = cntrl, 0% = full
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
region_name = 'grid'

# ROMS output location
outpath = '/data/project6/minnaho/opc_scenarios/bgc_flux/'

# roms var
var_name = 'npp' 
var_nc = 'var_int' 
cblabel = 'mmol m$^{-2}$ d$^{-1}$'

#year_month = 'Y1998_M04_06'
year_month1 = '1999'
year_month2 = '2016'

# scenario names 
exp = ['PNDN_only',
       'PNDN_only_realistic',
       'FNDN_only',
       'FNDN_only_realistic', 
       'pndn50', 
       'pndn50_realistic', 
       'pndn90', 
       'pndn90_realistic']

if year_month1 == 'fullts':
    title_exp = ['50% N\nReduction 98-99',
                 year_month2,
                 '85% N\nReduction 98-99',
                 year_month2,
                 '50% N\nReduction\n50% Recycle 98-99',
                 year_month2,
                 '50% N\nReduction\n90% Recycle 98-99',
                 year_month2]
    if year_month2 == 'fullts':
        title_exp = ['50% N\nReduction 98-99',
                     '16-17',
                     '85% N\nReduction 98-99',
                     '16-17',
                     '50% N\nReduction\n50% Recycle 98-99',
                     '16-17',
                     '50% N\nReduction\n90% Recycle 98-99',
                     '16-17']
else:
    title_exp = ['50% N\nReduction '+year_month1,
                 year_month2,
                 '85% N\nReduction '+year_month1,
                 year_month2,
                 '50% N\nReduction\n50% Recycle '+year_month1,
                 year_month2,
                 '50% N\nReduction\n90% Recycle '+year_month1,
                 year_month2]
    if year_month2 == 'fullts':
        title_exp = ['50% N\nReduction '+year_month1,
                     '16-17',
                     '85% N\nReduction '+year_month1,
                     '16-17',
                     '50% N\nReduction\n50% Recycle '+year_month1,
                     '16-17',
                     '50% N\nReduction\n90% Recycle '+year_month1,
                     '16-17']

#filest = 'int_avg_100m_50m_'
filest1 = 'avg_'+year_month1+'_int_avg_100m_50m_'
filest2 = 'avg_'+year_month2+'_int_avg_100m_50m_'

filest1_std = 'concat_'+year_month1+'_int_avg_100m_50m_'
filest2_std = 'concat_'+year_month2+'_int_avg_100m_50m_'

fileen = '_'+var_name+'.nc'


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

mask_mult[mask_mult==0] = np.nan

fpath = []
for e_i in range(len(exp)):
    if 'realistic' in exp[e_i]:
        fpath.append(outpath+filest2+exp[e_i]+fileen)
    else:
        fpath.append(outpath+filest1+exp[e_i]+fileen)

fpath_std = []
for e_i in range(len(exp)):
    if 'realistic' in exp[e_i]:
        fpath_std.append(outpath+filest2_std+exp[e_i]+fileen)
    else:
        fpath_std.append(outpath+filest1_std+exp[e_i]+fileen)
    

s2d = 86400

time_sh1 = Dataset(outpath+'concat_'+year_month1+'_int_avg_100m_50m_cntrl_initap'+fileen,'r').dimensions['time'].size
time_sh2 = Dataset(outpath+'concat_'+year_month2+'_int_avg_100m_50m_cntrl_2012_2017'+fileen,'r').dimensions['time'].size

# get cntrl and loads1617 to compare to
#cntrl_var = np.squeeze(Dataset(outpath+filest+'cntrl'+fileen,'r').variables[var_nc])*mask_mult
cntrl_var_old = np.squeeze(Dataset(outpath+'avg_'+year_month1+'_int_avg_100m_50m_cntrl_initap'+fileen,'r').variables[var_nc])*mask_mult
cntrl_mean_old = np.nanmean(cntrl_var_old)*s2d

cntrl_std_old_rd = np.squeeze(Dataset(outpath+'concat_'+year_month1+'_int_avg_100m_50m_cntrl_initap'+fileen,'r').variables[var_nc])*mask_mult
cntrl_std_old = np.nanstd(np.nanmean(cntrl_std_old_rd,axis=(1,2)))*s2d

cntrl_old_25 = np.nanpercentile(cntrl_std_old_rd,25)*s2d
cntrl_old_75 = np.nanpercentile(cntrl_std_old_rd,75)*s2d
cntrl_old_med = np.nanmedian(cntrl_std_old_rd)*s2d

cntrl_var_new = np.squeeze(Dataset(outpath+'avg_'+year_month2+'_int_avg_100m_50m_cntrl_2012_2017'+fileen,'r').variables[var_nc])*mask_mult
cntrl_mean_new = np.nanmean(cntrl_var_new)*s2d

cntrl_std_new_rd = np.squeeze(Dataset(outpath+'concat_'+year_month2+'_int_avg_100m_50m_cntrl_2012_2017'+fileen,'r').variables[var_nc])*mask_mult
cntrl_std_new = np.nanstd(np.nanmean(cntrl_std_new_rd,axis=(1,2)))*s2d

cntrl_new_25 = np.nanpercentile(cntrl_std_new_rd,25)*s2d
cntrl_new_75 = np.nanpercentile(cntrl_std_new_rd,75)*s2d
cntrl_new_med = np.nanmedian(cntrl_std_new_rd)*s2d


fulll_var_old = np.squeeze(Dataset(outpath+'avg_'+year_month1+'_int_avg_100m_50m_loads1617'+fileen,'r').variables[var_nc])*mask_mult
fulll_mean_old = np.nanmean(fulll_var_old)*s2d

fulll_std_old_rd = np.squeeze(Dataset(outpath+'concat_'+year_month1+'_int_avg_100m_50m_loads1617'+fileen,'r').variables[var_nc])*mask_mult
fulll_std_old = np.nanstd(np.nanmean(fulll_std_old_rd,axis=(1,2)))*s2d

fulll_old_25 = np.nanpercentile(fulll_std_old_rd,25)*s2d
fulll_old_75 = np.nanpercentile(fulll_std_old_rd,75)*s2d
fulll_old_med = np.nanmedian(fulll_std_old_rd)*s2d

fulll_var_new = np.squeeze(Dataset(outpath+'avg_'+year_month2+'_int_avg_100m_50m_fulll_2012_2017'+fileen,'r').variables[var_nc])*mask_mult
fulll_mean_new = np.nanmean(fulll_var_new)*s2d

fulll_std_new_rd = np.squeeze(Dataset(outpath+'concat_'+year_month2+'_int_avg_100m_50m_fulll_2012_2017'+fileen,'r').variables[var_nc])*mask_mult
fulll_std_new = np.nanstd(np.nanmean(fulll_std_new_rd,axis=(1,2)))*s2d

fulll_new_25 = np.nanpercentile(fulll_std_new_rd,25)*s2d
fulll_new_75 = np.nanpercentile(fulll_std_new_rd,75)*s2d
fulll_new_med = np.nanmedian(fulll_std_new_rd)*s2d

#l1617_var = np.squeeze(Dataset(outpath+filest+'loads1617'+fileen,'r').variables[var_nc])*mask_mult
#l1617_mean = np.nanmean(l1617_var)*s2d
#
#potwr_var = np.squeeze(Dataset(outpath+filest+'POTW'+fileen,'r').variables[var_nc])*mask_mult
#potwr_mean = np.nanmean(potwr_var)*s2d

# difference between l1617 and POTW only run = rivers 
#potwd_var = np.squeeze(Dataset(outpath+filest+'loads1617-POTW'+fileen,'r').variables[var_nc])*mask_mult
#potwd_mean = np.nanmean(potwd_var)*s2d

figw = 14
#figw = 10
figh = 4

axis_tick_size = 14
msize_s = 50
msize_l = 200

col1 = 'limegreen'
col2 = 'purple'

lw = 4

savename = var_name+'_'+year_month1+'_'+year_month2+'_'+region_name+'_'+exp[-1]+'_boxplot_cntrl_initap.png'
fig,ax = plt.subplots(1,1,figsize=[figw,figh])

roms_concat1 = np.ones((int(len(exp)/2),time_sh1))*np.nan
roms_concat2 = np.ones((int(len(exp)/2),time_sh2))*np.nan

i_1 = 0
i_2 = 0
for n_i in range(len(exp)):

    #average
    roms_var = np.squeeze(Dataset(fpath[n_i],'r').variables[var_nc])*mask_mult
    roms_neg_raw = np.nanmean(roms_var)*s2d

    # concatenated array
    roms_concat_rd = np.squeeze(Dataset(fpath_std[n_i],'r').variables[var_nc])*mask_mult
    roms_avg = np.nanmean(roms_concat_rd,axis=(1,2))*s2d

    #print(exp[n_i],'up to '+str(np.nanmin(((roms_var-l1617_var)/l1617_var)*100))+'% decrease in productivity')
    #print(exp[n_i],'up to '+str(np.nanmax(((roms_var-l1617_var)/l1617_var)*100))+'% increase in productivity')

    # calculate percent of full-cntrl
    #if 'real' in exp[n_i]:
    #    roms_neg = ((roms_neg_raw-fulll_mean_new)/(fulll_mean_new-cntrl_mean_new))*100
    #    roms_neg_std = ((roms_std-fulll_mean_new)/(fulll_mean_new-cntrl_mean_new))*100
    #    ax.bar(n_i,roms_neg,yerr=roms_neg_std,color='gray')
    #else:
    #    roms_neg = ((roms_neg_raw-fulll_mean_old)/(fulll_mean_old-cntrl_mean_old))*100
    #    roms_neg_std = ((roms_std-fulll_mean_old)/(fulll_mean_old-cntrl_mean_old))*100
    #    ax.bar(n_i,roms_neg,yerr=roms_neg_std,color='white',edgecolor='k')


    if 'real' in exp[n_i]:
        roms_concat2[i_2] = roms_avg
        ax.scatter(n_i,cntrl_new_med,marker='s',s=msize_s,color=col2)
        ax.scatter(n_i,cntrl_new_25,marker='_',s=msize_l,linewidth=lw,color=col2)
        ax.scatter(n_i,cntrl_new_75,marker='_',s=msize_l,linewidth=lw,color=col2)
        ax.scatter(n_i,fulll_new_med,marker='s',s=msize_s,color=col1)
        ax.scatter(n_i,fulll_new_25,marker='_',s=msize_l,linewidth=lw,color=col1)
        ax.scatter(n_i,fulll_new_75,marker='_',s=msize_l,linewidth=lw,color=col1)
        i_2 += 1
    else:
        roms_concat1[i_1] = roms_avg
        ax.scatter(n_i,cntrl_old_med,marker='s',s=msize_s,color=col2)
        ax.scatter(n_i,cntrl_old_25,marker='_',s=msize_l,linewidth=lw,color=col2)
        ax.scatter(n_i,cntrl_old_75,marker='_',s=msize_l,linewidth=lw,color=col2)
        ax.scatter(n_i,fulll_old_med,marker='s',s=msize_s,color=col1)
        ax.scatter(n_i,fulll_old_25,marker='_',s=msize_l,linewidth=lw,color=col1)
        ax.scatter(n_i,fulll_old_75,marker='_',s=msize_l,linewidth=lw,color=col1)
        i_1 += 1

    print(exp[n_i],'median',str(np.nanmedian(roms_avg)))
    print(exp[n_i],'avg',str(roms_neg_raw))

print('cntrl old median',str(cntrl_old_med))
print('cntrl old avg',str(cntrl_mean_old))
print('cntrl new median',str(cntrl_new_med))
print('cntrl new avg',str(cntrl_mean_new))
print('fulll old median',str(fulll_old_med))
print('fulll old avg',str(fulll_mean_old))
print('fulll new median',str(fulll_new_med))
print('fulll new avg',str(fulll_mean_new))

ax.boxplot(roms_concat1.T,positions=range(0,len(exp),2))
ax.boxplot(roms_concat2.T,positions=range(1,len(exp),2))

#ax.set_ylim(bottom=-100,top=10)
ax.set_xticks(range(len(exp)))

#ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*cntrl_mean_old,linestyle='-',color=col2)
#ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*cntrl_old_med,linestyle='-',color=col2)
#ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*cntrl_old_25,linestyle='--',color=col2)
#ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*cntrl_old_75,linestyle='--',color=col2)
#
##ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*cntrl_mean_new,linestyle='-',color='navy')
#ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*cntrl_new_med,linestyle='-',color='navy')
#ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*cntrl_new_25,linestyle='--',color='navy')
#ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*cntrl_new_75,linestyle='--',color='navy')
#
##ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*fulll_mean_old,linestyle='-',color=col1)
#ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*fulll_old_med,linestyle='-',color=col1)
#ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*fulll_old_25,linestyle='--',color=col1)
#ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*fulll_old_75,linestyle='--',color=col1)
#
##ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*fulll_mean_new,linestyle='-',color='orange')
#ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*fulll_new_med,linestyle='-',color='orange')
#ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*fulll_new_25,linestyle='--',color='orange')
#ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*fulll_new_75,linestyle='--',color='orange')


ax.set_xticklabels(title_exp,fontsize=12)
ax.set_xlabel('Scenario',fontsize=axis_tick_size)
ax.set_ylabel('NPP '+cblabel,fontsize=axis_tick_size)
ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)

fig.savefig(savepath+savename,bbox_inches='tight')
print(savename)
#plt.close()

