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

region_name = 'offshore'

outpath = './budget_ww/'
savepath = './figs/'

# scenario ANTH and CTRL
# put cntrl first in exp list
anth = 'L2SCB_AP'
cntrl = 'L2SCB'

exp = [cntrl,
       anth,
       'PNDN_only_realistic',
       'FNDN_only_realistic',
       'pndn50_realistic',
       'pndn90_realistic']

title_exp = ['',
             '',
             '50% N Reduction',
             '85% N Reduction',
             '50% N Reduction\n50% Recycle',
             '50% N Reduction\n90% Recycle']

varn = 'N'
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

mask_mult[mask_mult==0] = np.nan
s2d = 86400

figw = 16
figh = 6
axfont = 14

months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

# calculate BGC and uptake terms
clim_bgc = np.ones((len(exp),12))*np.nan
clim_bgc_std = np.ones((len(exp),12))*np.nan

clim_upt = np.ones((len(exp),12))*np.nan
clim_upt_std = np.ones((len(exp),12))*np.nan

clim_remin = np.ones((len(exp),12))*np.nan
clim_remin_std = np.ones((len(exp),12))*np.nan


for e_i in range(len(exp)):
    # read in file
    matr = h5py.File(outpath+fst+exp[e_i]+dep+'/'+matn+'.mat','r')
    data = matr.get(matn)

    # get dates
    datemat = np.squeeze(h5py.File(outpath+fst+exp[e_i]+dep+'/'+matc+'.mat','r').get(matc)['date'])

    # get variables
    nitrif = np.squeeze(data['NITRIF'])
    denitr = np.squeeze(data['DENIT'])
    seddenitr = np.squeeze(data['SED_DENITR'])
    no3up = np.squeeze(data['PHOTO_NO3'])
    nh4up = np.squeeze(data['PHOTO_NH4'])
    donre = np.squeeze(data['DON_REMIN'])
    pocre = np.squeeze(data['POC_REMIN'])
    sedre = np.squeeze(data['SED_REMIN'])
    biore = np.squeeze(data['BIOLOGICAL_RELEASE'])
    ammox = np.squeeze(data['AMMOX'])

    # calculate terms
    bgc_no3_calc = nitrif-denitr-seddenitr-no3up
    bgc_nh4_calc = donre-nh4up+pocre+biore+sedre-ammox
    dinuptake = no3up+nh4up

    remin_no3 = nitrif-denitr-seddenitr
    remin_nh4 = donre+pocre+biore+sedre-ammox
    
    # average over mask
    if exp[e_i] == cntrl:
        bgc_cntrl = np.nanmean(((bgc_no3_calc+bgc_nh4_calc)*mask_mult),axis=(1,2))*s2d
        upt_cntrl = np.nanmean(dinuptake*mask_mult,axis=(1,2))*s2d
        remin_cntrl = np.nanmean(((remin_no3+remin_nh4)*mask_mult),axis=(1,2))*s2d
    
    else:
        bgc = (np.nanmean(((bgc_no3_calc+bgc_nh4_calc)*mask_mult),axis=(1,2))*s2d)-bgc_cntrl
        upt = (np.nanmean(dinuptake*mask_mult,axis=(1,2))*s2d)-upt_cntrl
        remin = (np.nanmean(((remin_no3+remin_nh4)*mask_mult),axis=(1,2))*s2d)-remin_cntrl
        

        # convert matlab time
        dt = pd.to_datetime(datemat-719529, unit='D')
        for m_i in range(1,13):
            clim_bgc[e_i,m_i-1] = np.nanmean(bgc[np.where(dt.month==m_i)[0]])
            clim_bgc_std[e_i,m_i-1] = np.nanstd(bgc[np.where(dt.month==m_i)[0]])

            clim_upt[e_i,m_i-1] = np.nanmean(upt[np.where(dt.month==m_i)[0]])
            clim_upt_std[e_i,m_i-1] = np.nanstd(upt[np.where(dt.month==m_i)[0]])

            clim_remin[e_i,m_i-1] = np.nanmean(remin[np.where(dt.month==m_i)[0]])
            clim_remin_std[e_i,m_i-1] = np.nanstd(remin[np.where(dt.month==m_i)[0]])

        # plot
        if exp[e_i] != anth:
            fig,ax = plt.subplots(3,1,figsize=[figw,figh])

            # plot scenarios as bars
            ax[0].bar(range(1,13),clim_bgc[e_i],yerr=clim_bgc_std[e_i],color='white',edgecolor='k',capsize=4)
            ax[0].set_title('Remineralization - Uptake',fontsize=axfont)

            ax[1].bar(range(1,13),clim_upt[e_i],yerr=clim_upt_std[e_i],color='green',edgecolor='k',capsize=4)
            ax[1].set_title('Uptake',fontsize=axfont)

            ax[2].bar(range(1,13),clim_remin[e_i],yerr=clim_remin_std[e_i],color='blue',edgecolor='k',capsize=4)
            ax[2].set_title('Remineralization',fontsize=axfont)

            # plot anth as line
            ax[0].plot(range(1,13),clim_bgc[exp.index(anth)])
            ax[0].fill_between(range(1,13),clim_bgc[exp.index(anth)]+clim_bgc_std[exp.index(anth)],clim_bgc[exp.index(anth)]-clim_bgc_std[exp.index(anth)],alpha=0.3)

            ax[1].plot(range(1,13),clim_upt[exp.index(anth)])
            ax[1].fill_between(range(1,13),clim_upt[exp.index(anth)]+clim_upt_std[exp.index(anth)],clim_upt[exp.index(anth)]-clim_upt_std[exp.index(anth)],alpha=0.3)

            ax[2].plot(range(1,13),clim_remin[exp.index(anth)])
            ax[2].fill_between(range(1,13),clim_remin[exp.index(anth)]+clim_remin_std[exp.index(anth)],clim_remin[exp.index(anth)]-clim_remin_std[exp.index(anth)],alpha=0.3)

            for ax_i in range(ax.shape[0]):
                ax[ax_i].set_ylabel('mmol N m$^{-3}$ d$^{-1}$',fontsize=axfont)
                ax[ax_i].tick_params(axis='both',which='major',labelsize=axfont)
                ax[ax_i].set_xticks(range(1,13))


            ax[0].set_xticklabels('')
            ax[1].set_xticklabels('')
            ax[2].set_xticklabels(months,fontsize=axfont)
            if region_name == 'grid':
                ax[0].set_ylim(bottom=-.75,top=.5)
                ax[1].set_ylim(bottom=-.5,top=1.75)
                ax[2].set_ylim(bottom=-.5,top=.75)
            if region_name == 'coast':
                ax[0].set_ylim(bottom=-3,top=.5)
                ax[1].set_ylim(bottom=-.5,top=5)
                ax[2].set_ylim(bottom=-.5,top=2.5)
            if region_name == 'offshore':
                ax[0].set_ylim(bottom=-.5,top=.5)
                ax[1].set_ylim(bottom=-.5,top=1.5)
                ax[2].set_ylim(bottom=-.5,top=1)
            #fig.suptitle(title_exp[e_i],fontsize=axfont)
            savename = 'massb_clim_'+varn+dep+'_'+exp[e_i]+'_'+region_name+'.png'
            fig.savefig(savepath+savename,bbox_inches='tight')

