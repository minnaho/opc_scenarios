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

region_name = 'grid'

outpath = './budget_ww/'
savepath = './figs/'

timep = 'julnov'

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
    #jul-nov
    if timep == 'julnov':
        timest1 = '2016-07-01'
        timeen1 = '2016-12-01'
        timest2 = '2017-07-01'
        timeen2 = '2017-12-01'

        timest3 = '1998-07-01'
        timeen3 = '1998-11-30'
        timest4 = '1999-07-01'
        timeen4 = '1999-11-30'

    #full timeperiod
    if timep == 'fullts':
        timest1 = '2015-11-01'
        timeen1 = '2016-12-01'
        timest2 = '2016-11-01'
        timeen2 = '2017-12-01'

        timest3 = '1997-11-01'
        timeen3 = '1998-11-30'
        timest4 = '1998-11-01'
        timeen4 = '1999-11-30'

    #sep-nov
    if timep == 'sepnov':
        timest1 = '2016-09-01'
        timeen1 = '2016-12-01'
        timest2 = '2017-09-01'
        timeen2 = '2017-12-01'

        timest3 = '1998-09-01'
        timeen3 = '1998-12-01'
        timest4 = '1998-09-01'
        timeen4 = '1999-12-01'

# scenario ANTH and CTRL
anth1 = 'L2SCB_AP'
cntrl1 = 'L2SCB'
anth2 = 'loads1617'
cntrl2 = 'cntrl_initap'

if sce == 'recy':
    exp1 = [cntrl1,
           anth1,
           'PNDN_only_realistic',
           'pndn50_realistic',
           'pndn90_realistic',
           'FNDN_only_realistic'
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
    title_exp1 = [cntrl1,
                 anth1,
                 '50% N Red.',
                 '50% N Red.\n50% Recy.',
                 '50% N Red.\n90% Recy.',
                 '85% N Red.',
                 '85% N Red.\n50% Recy.',
                 '85% N Red.\n90% Recy.']
    title_exp2=  [
                 cntrl2,
                 anth2,
                 '50% N Red.',
                 '50% N Red.\n50% Recy.',
                 '50% N Red.\n90% Recy.',
                 '85% N Red.',
                ]

if sce == 'nman':
    exp1 = [cntrl1,
           anth1,
           'PNDN_only_realistic',
           'FNDN_only_realistic'
                                ]
    exp2 = [cntrl2,
           anth2,
           'PNDN_only',
           'FNDN_only'
                            ]
    title_exp = [
                 '50% N Red.',
                 '85% N Red.',
                 '50% N Red.',
                 '85% N Red.'
                ]


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
    regtitle = 'Offshore'
    mask_temp = np.copy(mask_cst)
    mask_temp[mask_temp==1] = 2
    mask_temp[np.isnan(mask_temp)] = 1
    mask_temp[mask_temp==2] = 0
    mask_temp = mask_temp*mask_nc
    mask_temp[:,:20] = np.nan
    mask_temp[:20,:] = np.nan
    mask_temp[-20:,:] = np.nan
    mask_mult = mask_temp
    #regtitle = 'Offshore'
    #mask7[mask9==1] = 1
    #mask_mult = mask7


mask_mult[mask_mult==0] = np.nan
s2d = 86400

if sce == 'recy':
    figw = 16
if sce == 'nman':
    figw = 8

figh = 4
axfont = 14

bmin = 20
bmax = 200

nbins = 500
bins_p = np.linspace(bmin,bmax,nbins+1)

# pdf function
n_p1 = np.zeros((len(exp1),nbins+1))
n_p2 = np.zeros((len(exp2),nbins+1))

flt1 = np.zeros((len(exp2)))
flt2 = np.zeros((len(exp2)))

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
    respir_calc = dinuptake*(117./16)*s2d # too lazy to change var names below
    respir_calc[respir_calc==0] = np.nan

    # average over mask
    #if exp1[e_i] == cntrl1:
    #    respir_cntrl1 = np.copy(respir_calc)*s2d
    #    respir_cntrl1[respir_cntrl1==0] = np.nan
    #    n,bins,patch = plt.hist(respir_cntrl1.flatten(),bins=nbins,range=([bmin,bmax]))

    #else:
    n_e,bins,patch = plt.hist(respir_calc.flatten(),bins=nbins,range=([bmin,bmax]))
    n_e = np.append(n_e,0)
    n_p1[e_i] = n_e

    flt1[e_i] = np.where(~np.isnan(respir_calc.flatten()))[0].shape[0]
    
        #respir_avg1 = (respir_calc*s2d)-respir_cntrl1

        #respir_avg1[respir_avg1==0] = np.nan

        ## convert matlab time
        #dt = pd.to_datetime(datemat-719529, unit='D')

        ## calculate avg and std over each year
        #clim_res1[e_i,0] = np.nanmean(respir_avg1[((dt>timest1)&(dt<timeen1))])
        #clim_res1[e_i,1] = np.nanmean(respir_avg1[((dt>timest2)&(dt<timeen2))])
        #if exp1[e_i] == anth1:
        #    print(exp1[e_i],str(clim_res1[e_i,0]))

for e_i in range(len(exp2)):
    # read in file
    matr = h5py.File(outpath+fst+exp2[e_i]+dep+'/'+matn+'.mat','r')
    data = matr.get(matn)

    # get dates
    datemat = np.squeeze(h5py.File(outpath+fst+exp2[e_i]+dep+'/'+matc+'.mat','r').get(matc)['date'])

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
    respir_calc = dinuptake*(117./16)*s2d # too lazy to change var names below
    respir_calc[respir_calc==0] = np.nan
    n_e,bins,patch = plt.hist(respir_calc.flatten(),bins=nbins,range=([bmin,bmax]))
    n_e = np.append(n_e,0)
    n_p2[e_i] = n_e
    flt2[e_i] = np.where(~np.isnan(respir_calc.flatten()))[0].shape[0]
    
    # average over mask
    #if exp2[e_i] == cntrl2:
    #    respir_cntrl2 = np.nanmean(((respir_calc)*mask_mult),axis=(1,2))*s2d
    #    respir_cntrl2[respir_cntrl2==0] = np.nan

    #else:
    #
    #    respir_avg2 = (np.nanmean(((respir_calc)*mask_mult),axis=(1,2))*s2d)-respir_cntrl2

    #    respir_avg2[respir_avg2==0] = np.nan

    #    # convert matlab time
    #    dt = pd.to_datetime(datemat-719529, unit='D')

    #    # calculate avg and std over each year
    #    clim_res2[e_i,0] = np.nanmean(respir_avg2[((dt>timest3)&(dt<timeen3))])
    #    clim_res2[e_i,1] = np.nanmean(respir_avg2[((dt>timest4)&(dt<timeen4))])

if sce == 'recy':
    ind1 = list(range(6,10))
    ind2 = list(range(0,6))
if sce == 'nman':
    ind2 = list(range(0,2))
    ind1 = list(range(2,4))

# plot
fig,ax = plt.subplots(1,1,figsize=[figw,figh])

for p_i in range(len(exp1)):
    ax.plot(bins_p,n_p1[p_i]/flt1[p_i],label=title_exp1[e_i])
for p_i in range(len(exp2)):
    ax.plot(bins_p,n_p2[p_i]/flt2[p_i],label=title_exp2[e_i])

savename = 'pdf_uptakeC_'+region_name+'_'+timep+'_'+sce+'.png'
fig.savefig(savepath+savename,bbox_inches='tight')

