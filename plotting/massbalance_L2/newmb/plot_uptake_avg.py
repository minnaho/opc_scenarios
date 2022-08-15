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

sce = 'recy'


if region_name == 'coast':
    timest1 = '2015-11-01'
    timeen1 = '2016-11-01'
    timest2 = '2016-11-01'
    timeen2 = '2017-11-01'

    timest3 = '1997-11-01'
    timeen3 = '1998-11-01'
    timest4 = '1998-11-01'
    timeen4 = '1999-11-01'

if region_name == 'offshore' or region_name == 'grid':
    timest1 = '2016-08-01'
    timeen1 = '2016-12-31'
    timest2 = '2017-08-01'
    timeen2 = '2017-12-31'

    timest3 = '1998-08-01'
    timeen3 = '1998-12-31'
    timest4 = '1999-08-01'
    timeen4 = '1999-12-31'

# scenario ANTH and CTRL
anth2 = 'loads1617'
cntrl2 = 'cntrl_initap'
anth1 = 'L2SCB_AP'
cntrl1 = 'L2SCB'

exp1 = [cntrl1,
       anth1,
       'PNDN_only_realistic',
       'FNDN_only_realistic',
       'pndn50_realistic',
       'pndn90_realistic'
                            ]

exp2 = [cntrl2,
       anth2,
       'PNDN_only',
       'pndn50',
       'pndn90',
       'FNDN_only',
       'fndn50',
       'fndn90'
                            ]



#exp = [cntrl1,
#       anth1,
#       cntrl2,
#       anth2,
#       'PNDN_only_realistic',
#       'FNDN_only_realistic',
#       'pndn50_realistic',
#       'pndn90_realistic'
#       'PNDN_only',
#       'FNDN_only',
#       'pndn50',
#       'pndn90'
#       'fndn50',
#       'fndn90'
#                            ]

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
    #mask_temp = np.copy(mask_cst)
    #mask_temp[mask_temp==1] = 2
    #mask_temp[np.isnan(mask_temp)] = 1
    #mask_temp[mask_temp==2] = 0
    #mask_temp = mask_temp*mask_nc
    #mask_temp[:,:20] = np.nan
    #mask_temp[:20,:] = np.nan
    #mask_temp[-20:,:] = np.nan
    #mask_mult = mask_temp
    regtitle = 'Offshore'
    mask7[mask9==1] = 1
    mask_mult = mask7
    regtitle = 'Offshore'


mask_mult[mask_mult==0] = np.nan
s2d = 86400

figw = 16
figh = 4
axfont = 14

months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

# calculate BGC and uptake terms
clim_res1 = np.ones((len(exp1),2))*np.nan

for e_i in range(len(exp1)):
    # read in file
    matr = h5py.File(outpath+fst+exp1[e_i]+dep+'/'+matn+'.mat','r')
    data = matr.get(matn)

    # get dates
    datemat = np.squeeze(h5py.File(outpath+fst+exp1[e_i]+dep+'/'+matc+'.mat','r').get(matc)['date'])

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

    dinuptake = no3up+nh4up
    respir_calc = dinuptake # too lazy to change var names below

    # average over mask
    if exp1[e_i] == cntrl1:
        respir_cntrl1 = np.nanmean(((respir_calc)*mask_mult),axis=(1,2))*s2d
        respir_cntrl1[respir_cntrl1==0] = np.nan

    else:
    
        respir_avg1 = (np.nanmean(((respir_calc)*mask_mult),axis=(1,2))*s2d)-respir_cntrl1

        respir_avg1[respir_avg1==0] = np.nan

        # convert matlab time
        dt = pd.to_datetime(datemat-719529, unit='D')

        # calculate avg and std over each year
        clim_res1[e_i,0] = np.nanmean(respir_avg1[((dt>timest1)&(dt<timeen1))])
        clim_res1[e_i,1] = np.nanmean(respir_avg1[((dt>timest2)&(dt<timeen2))])
        if exp1[e_i] == anth1:
            print(exp1[e_i],str(clim_res1[e_i,0]))

clim_res_std1 = np.nanstd(clim_res1,axis=1)[2:]
plt1 = np.nanmean(clim_res1,axis=1)[2:]
anthmean1 = np.nanmean(clim_res1[1])
anthstd1 = np.nanstd(clim_res1[1])


clim_res2 = np.ones((len(exp2),2))*np.nan

for e_i in range(len(exp2)):
    # read in file
    matr = h5py.File(outpath+fst+exp2[e_i]+dep+'/'+matn+'.mat','r')
    data = matr.get(matn)

    # get dates
    datemat = np.squeeze(h5py.File(outpath+fst+exp2[e_i]+dep+'/'+matc+'.mat','r').get(matc)['date'])

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
    if exp2[e_i] == cntrl2:
        respir_cntrl2 = np.nanmean(((respir_calc)*mask_mult),axis=(1,2))*s2d
        respir_cntrl2[respir_cntrl2==0] = np.nan

    else:
    
        respir_avg2 = (np.nanmean(((respir_calc)*mask_mult),axis=(1,2))*s2d)-respir_cntrl2

        respir_avg2[respir_avg2==0] = np.nan

        # convert matlab time
        dt = pd.to_datetime(datemat-719529, unit='D')

        # calculate avg and std over each year
        clim_res2[e_i,0] = np.nanmean(respir_avg2[((dt>timest3)&(dt<timeen3))])
        clim_res2[e_i,1] = np.nanmean(respir_avg2[((dt>timest4)&(dt<timeen4))])

clim_res_std2 = np.nanstd(clim_res2,axis=1)[2:]
plt2 = np.nanmean(clim_res2,axis=1)[2:]

anthmean2 = np.nanmean(clim_res2[1])
anthstd2 = np.nanstd(clim_res2[1])

exit()

if sce == 'recy':
    ind1 = list(range(0,6))
    ind2 = list(range(6,10))
if sce == 'nman':
    ind1 = list(range(0,2))
    ind2 = list(range(2,4))

# plot
fig,ax = plt.subplots(1,1,figsize=[figw,figh])

# plot scenarios as bars
ax.bar(ind1,plt1,yerr=clim_res_std1,color='white',edgecolor='k',capsize=4)
ax.bar(ind2,plt2,yerr=clim_res_std2,color='gray',capsize=4)
ax.set_title(regtitle,fontsize=axfont)

# plot anth as line
ax.plot(ind1,np.ones((len(ind1)))*anthmean1)
ax.fill_between(ind1,np.ones((len(ind1)))*anthmean1+anthstd1,np.ones((len(ind1)))*anthmean1+anthstd1,alpha=0.3)

ax.plot(ind2,np.ones((len(ind2)))*anthmean2)
ax.fill_between(ind2,np.ones((len(ind2)))*anthmean2+anthstd2,np.ones((len(ind2)))*anthmean2+anthstd2,alpha=0.3)

ax.set_ylabel('Uptake\nmmol N m$^{-3}$ d$^{-1}$',fontsize=axfont)
ax.tick_params(axis='both',which='major',labelsize=axfont)

#ax.set_xticks(range(1,13))
#ax.set_xticklabels(months,fontsize=axfont)

#if region_name == 'grid':
#    ax.set_ylim(bottom=-.75,top=.5)
#if region_name == 'coast':
#    ax.set_ylim(bottom=-3,top=.5)
#if region_name == 'offshore':
#    ax.set_ylim(bottom=-.5,top=.5)

savename = 'massb_avg_'+varn+dep+'_uptake_'+region_name+'.png'
fig.savefig(savepath+savename,bbox_inches='tight')

