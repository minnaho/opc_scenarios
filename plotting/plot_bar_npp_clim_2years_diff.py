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
region_name = 'grid'

# ROMS output location
outpath = '/data/project6/minnaho/opc_scenarios/bgc_flux/'

# roms var
var_name = 'npp' 
var_nc = 'var_int' 
cblabel = 'mmol m$^{-2}$ d$^{-1}$'

filest = 'int_avg_100m_50m_'

#year_month = 'Y1998_M04_06'
year_month1 = 'fullts'
year_month2 = ''

# scenario names 
exp = ['PNDN_only',
       'PNDN_only_realistic',
       'FNDN_only',
       'FNDN_only_realistic', 
       'pndn50', 
       'pndn50_realistic', 
       'pndn90', 
       'pndn90_realistic',
       'fndn50',
       'fndn90']

# mask
mask_nc = l2grid.mask_nc


if year_month1 == 'fullts':
    title_exp = ['50% N\nReduction\n98-99',
                 year_month2,
                 '85% N\nReduction\n98-99',
                 year_month2,
                 '50% N\nReduction\n50% Recycle\n98-99',
                 year_month2,
                 '50% N\nReduction\n90% Recycle\n98-99',
                 year_month2,
                 '85% N\nReduction\n50% Recycle\n98-99',
                 '85% N\nReduction\n90% Recycle\n98-99']
    if year_month2 == 'fullts':
        title_exp = ['50% N\nReduction\n98-99',
                     '16-17',
                     '85% N\nReduction\n98-99',
                     '16-17',
                     '50% N\nReduction\n50% Recycle\n98-99',
                     '16-17',
                     '50% N\nReduction\n90% Recycle\n98-99',
                     '16-17',
                     '85% N\nReduction\n50% Recycle\n98-99',
                     '85% N\nReduction\n90% Recycle\n98-99']
else:
    title_exp = ['50% N\nReduction\n'+year_month1,
                 year_month2,
                 '85% N\nReduction\n'+year_month1,
                 year_month2,
                 '50% N\nReduction\n50% Recycle\n'+year_month1,
                 year_month2,
                 '50% N\nReduction\n90% Recycle\n'+year_month1,
                 year_month2,
                 '85% N\nReduction\n50% Recycle\n'+year_month1,
                 '85% N\nReduction\n90% Recycle\n'+year_month1]
    if year_month2 == 'fullts':
        title_exp = ['50% N\nReduction\n'+year_month1,
                     '16-17',
                     '85% N\nReduction\n'+year_month1,
                     '16-17',
                     '50% N\nReduction\n50% Recycle\n'+year_month1,
                     '16-17',
                     '50% N\nReduction\n90% Recycle\n'+year_month1,
                     '16-17',
                     '85% N\nReduction\n50% Recycle\n'+year_month1,
                     '85% N\nReduction\n90% Recycle\n'+year_month1]


filest1 = 'avg_'+year_month1+'_int_avg_100m_50m_'
filest2 = 'avg_'+year_month2+'_int_avg_100m_50m_'
fileen = '_'+var_name+'.nc'

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

# seconds to days
s2d = 86400

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

#c_i = 0
#t_i = 0
clim = np.ones((len(exp),12))*np.nan
cstd = np.ones((len(exp),12))*np.nan
# per month
for m_i in range(1,13):
    for e_i in range(len(exp)):
        #t_i += 1
        flist = glob.glob(outpath+filest+exp[e_i]+'_'+var_name+'*M%02d'%m_i+'.nc')
        #print('len flist',str(len(flist)))
        if len(flist)>1:
            #c_i += 1
            temparr = np.ones((len(flist),mask_nc.shape[0],mask_nc.shape[1]))*np.nan
            for f_i in range(len(flist)):
                temparr[f_i] = (np.squeeze(Dataset(flist[f_i],'r').variables[var_nc])*mask_mult)
                temparr[temparr==0] = np.nan
                if 'real' in exp[e_i]:
                    cntrlm = np.squeeze(Dataset(outpath+filest+'cntrl_2012_2017_'+var_name+'_'+flist[f_i][flist[f_i].index('Y'):],'r').variables[var_nc])
                else:
                    cntrlm = np.squeeze(Dataset(outpath+filest+'cntrl_initap_'+var_name+'_'+flist[f_i][flist[f_i].index('Y'):],'r').variables[var_nc])
                cntrlm[cntrlm==0] = np.nan
                tempdiff = temparr - cntrlm
            cstd[e_i,m_i-1] = np.nanstd(tempdiff)
            clim[e_i,m_i-1] = np.nanmean(tempdiff)
#        else:
#            tempnc = np.squeeze(Dataset(flist[0],'r').variables[var_nc])*mask_mult
#            tempnc[tempnc==0] = np.nan
#            if 'real' in exp[e_i]:
#                cntrlm = np.squeeze(Dataset(outpath+filest+'cntrl_2012_2017_'+var_name+'_'+flist[0][flist[0].index('Y'):],'r').variables[var_nc])
#            else:
#                cntrlm = np.squeeze(Dataset(outpath+filest+'cntrl_initap_'+var_name+'_'+flist[0][flist[0].index('Y'):],'r').variables[var_nc])
#            tempdiff = tempnc - cntrlm
#            clim[e_i,m_i-1] = np.nanmean(tempdiff)

clim = clim*s2d
cstd = cstd*s2d

            
                
                

fpath = []
for e_i in range(len(exp)):
    if 'realistic' in exp[e_i]:
        fpath.append(outpath+filest2+exp[e_i]+fileen)
    else:
        fpath.append(outpath+filest1+exp[e_i]+fileen)



#figw = 14
figw = 16
figh = 6

axis_tick_size = 14

savename = var_name+'_'+year_month1+'_'+year_month2+'_'+region_name+'_'+exp[-1]+'_bar_norm_compoldnew_recy_cntrlinitap_diff.png'
fig,ax = plt.subplots(1,1,figsize=[figw,figh])

for n_i in range(len(exp)):

    roms_var = np.squeeze(Dataset(fpath[n_i],'r').variables[var_nc])*mask_mult
    roms_var[roms_var==0] = np.nan
    roms_neg = np.nanmean(roms_var)*s2d

    print(exp[n_i],str(roms_neg))
    #print(exp[n_i],'up to '+str(np.nanmin(((roms_var-l1617_var)/l1617_var)*100))+'% decrease in productivity')
    #print(exp[n_i],'up to '+str(np.nanmax(((roms_var-l1617_var)/l1617_var)*100))+'% increase in productivity')
    if 'real' in exp[n_i]:
        ax.bar(n_i,roms_neg-cntrl_mean_new,color='gray')
    else:
        ax.bar(n_i,roms_neg-cntrl_mean_old,color='white',edgecolor='k')

#ax.set_ylim(bottom=0,top=(fulll_mean_new-cntrl_mean_new)+0.5)
#ax.set_ylim(bottom=0)
ax.set_xticks(range(len(exp)))
#ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*cntrl_mean_old,linestyle='--',color='purple')
#ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*cntrl_mean_new,linestyle=':',color='purple')
ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*(fulll_mean_old-cntrl_mean_old),linestyle='--',color='green')
ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*(fulll_mean_new-cntrl_mean_new),linestyle=':',color='green')
#ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*l1617_mean,linestyle='--',color='green')
#ax.plot(range(-1,len(exp)+1),np.ones((len(exp)+2))*potwd_mean,linestyle='--',color='blue')

ax.set_xticklabels(title_exp,fontsize=12,rotation=90)
ax.set_xlabel('Scenario',fontsize=axis_tick_size)
ax.set_ylabel('Integrated NPP 100 m '+cblabel,fontsize=axis_tick_size)
ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)

fig.savefig(savepath+savename,bbox_inches='tight')
print(savename)
#plt.close()

