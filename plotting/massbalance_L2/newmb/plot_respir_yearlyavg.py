import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import l2grid as l2grid
import numpy as np
import h5py
from netCDF4 import Dataset
import pandas as pd
import matplotlib.pyplot as plt

plt.ion()

region_name = 'coast'

outpath = './budget_ww/'
savepath = './figs/'

# scenario ANTH and CTRL
# put cntrl first in exp list
anth1 = 'loads1617'
cntrl1 = 'cntrl_initap'
anth2 = 'L2SCB_AP'
cntrl2 = 'L2SCB'

exp = [cntrl1,
       anth1,
       cntrl2,
       anth2,
       'PNDN_only_realistic',
       'FNDN_only_realistic',
       'pndn50_realistic',
       'pndn90_realistic'
       'PNDN_only',
       'FNDN_only',
       'pndn50',
       'pndn90'
       'fndn50',
       'fndn90'
                            ]

title_exp = ['',
             '',
             '50% N Reduction',
             '85% N Reduction',
             '50% N Reduction\n50% Recycle',
             '50% N Reduction\n90% Recycle',
             '85% N Reduction\n50% Recycle',
             '85% N Reduction\n90% Recycle']

varn = 'O2'
matn = 'MATBGCF'
matc = 'MATVARC'

fst = 'outputs_'+varn+'_'
dep = '_0-200m'


mask_nc = l2grid.mask_nc

region_mask = Dataset('/data/project1/minnaho/make_masks/mask_scb.nc','r')
mask_cst = np.array(region_mask.variables['mask_coast'])
mask_cst[mask_cst==0] = np.nan

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
    regtitle = 'Offshore'
    '''
    # faycal's offshore mask7+mask9
    mask7[mask9==1] = 1
    #mask7[mask8==1] = 1 
    mask_mult = mask7


mask_mult[mask_mult==0] = np.nan
s2d = 86400

figw = 16
figh = 4
axfont = 14


# calculate BGC and uptake terms
clim_res = np.ones((len(exp),12))*np.nan
clim_res_std = np.ones((len(exp),12))*np.nan

for e_i in range(len(exp)):
    # read in file
    matr = h5py.File(outpath+fst+exp[e_i]+dep+'/'+matn+'.mat','r')
    data = matr.get(matn)

    # get dates
    datemat = np.squeeze(h5py.File(outpath+fst+exp[e_i]+dep+'/'+matc+'.mat','r').get(matc)['date'])

    # get variables
    loss = np.squeeze(data['LOSS'])
    graze = np.squeeze(data['GRAZE'])
    remin = np.squeeze(data['REMIN'])
    sedre = np.squeeze(data['SED_REMIN'])
    ammox = np.squeeze(data['AMMOX'])
    nit = np.squeeze(data['NIT'])

    # calculate terms
    respir_calc = loss+graze+remin+sedre+ammox+nit
    
    # average over mask
    if exp[e_i] == cntrl1:
        respir_cntrl1 = np.nanmean(((respir_calc)*mask_mult),axis=(1,2))*s2d
        respir_cntrl1[respir_cntrl1==0] = np.nan

    elif exp[e_i] == cntrl2:
        respir_cntrl2 = np.nanmean(((respir_calc)*mask_mult),axis=(1,2))*s2d
        respir_cntrl2[respir_cntrl2==0] = np.nan

    else:
    
        if 'realistic' in exp[e_i] :
            respir_avg = (np.nanmean(((respir_calc)*mask_mult),axis=(1,2))*s2d)-respir_cntrl2

        else:
            respir_avg = (np.nanmean(((respir_calc)*mask_mult),axis=(1,2))*s2d)-respir_cntrl1

        respir_avg[respir_avg==0] = np.nan

        # convert matlab time
        dt = pd.to_datetime(datemat-719529, unit='D')

        # calculate avg and std each year
        for m_i in range(1,13):
            clim_res[e_i,m_i-1] = np.nanmean(respir_avg[np.where(dt.month==m_i)[0]])
            clim_res_std[e_i,m_i-1] = np.nanstd(respir_avg[np.where(dt.month==m_i)[0]])

        # plot
        if exp[e_i] != anth:
            fig,ax = plt.subplots(1,1,figsize=[figw,figh])

            # plot scenarios as bars
            ax.bar(range(1,13),clim_res[e_i],yerr=clim_res_std[e_i],color='white',edgecolor='k',capsize=4)
            ax.set_title(title_exp[e_i],fontsize=axfont)

            # plot anth as line
            ax.plot(range(1,13),clim_res[exp.index(anth)])
            ax.fill_between(range(1,13),clim_res[exp.index(anth)]+clim_res_std[exp.index(anth)],clim_res[exp.index(anth)]-clim_res_std[exp.index(anth)],alpha=0.3)

            ax.set_ylabel('mmol O m$^{-3}$ d$^{-1}$',fontsize=axfont)
            ax.tick_params(axis='both',which='major',labelsize=axfont)
            ax.set_xticks(range(1,13))


            ax.set_xticklabels(months,fontsize=axfont)
            #if region_name == 'grid':
            #    ax.set_ylim(bottom=-.75,top=.5)
            #if region_name == 'coast':
            #    ax.set_ylim(bottom=-3,top=.5)
            #if region_name == 'offshore':
            #    ax.set_ylim(bottom=-.5,top=.5)
            #fig.suptitle(title_exp[e_i],fontsize=axfont)
            savename = 'massb_clim_'+varn+dep+'_'+exp[e_i]+'_'+region_name+'.png'
            fig.savefig(savepath+savename,bbox_inches='tight')

