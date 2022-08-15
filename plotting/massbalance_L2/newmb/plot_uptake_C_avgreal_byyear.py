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
    title_exp1 = [
                 '50% N Red.',
                 '50% N Red.\n50% Recy.',
                 '50% N Red.\n90% Recy.',
                 '85% N Red.'
                 ]
    title_exp2 = [
                 '50% N Red.',
                 '50% N Red.\n50% Recy.',
                 '50% N Red.\n90% Recy.',
                 '85% N Red.',
                 '85% N Red.\n50% Recy.',
                 '85% N Red.\n90% Recy.'
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
    title_exp1 = [
                 '50% N Red.',
                 '85% N Red.',
                 '50% N Red.',
                 '85% N Red.'
                ]
    title_exp2 = [
                 '50% N Red.',
                 '85% N Red.'
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


varn = 'N'
matn = 'MATBGCF'
matc = 'MATVARC'

fst = 'outputs_'+varn+'_'
dep = '_0-200m'


mask_nc = l2grid.mask_nc

region_mask = Dataset('/data/project1/minnaho/make_masks/mask_scb.nc','r')
mask_cst = np.array(region_mask.variables['mask_coast'])
mask_cst[mask_cst==0] = np.nan

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
    #mask_temp = np.copy(mask_cst)
    #mask_temp[mask_temp==1] = 2
    #mask_temp[np.isnan(mask_temp)] = 1
    #mask_temp[mask_temp==2] = 0
    #mask_temp = mask_temp*mask_nc
    #mask_temp[:,:20] = np.nan
    #mask_temp[:20,:] = np.nan
    #mask_temp[-20:,:] = np.nan
    #mask_mult = mask_temp
    #regtitle = 'Offshore'
    mask7[mask9==1] = 1
    mask_mult = mask7


mask_mult[mask_mult==0] = np.nan
s2d = 86400

if sce == 'recy':
    figw2 = 12
    figw1 = 8
if sce == 'nman':
    figw1 = 4
    figw2 = 4

figh = 4
axfont = 14

# calculate BGC and uptake terms
clim_res1 = np.ones((len(exp1),2))*np.nan
clim_std1 = np.ones((len(exp1),2))*np.nan

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
    respir_calc = dinuptake*(117./16) # too lazy to change var names below

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
        clim_std1[e_i,0] = np.nanstd(respir_avg1[((dt>timest1)&(dt<timeen1))])
        clim_std1[e_i,1] = np.nanstd(respir_avg1[((dt>timest2)&(dt<timeen2))])
        if exp1[e_i] == anth1:
            print(exp1[e_i],str(clim_res1[e_i,0]))

plt1 = np.nanmean(clim_res1,axis=1)[2:]
anthmean1 = np.nanmean(clim_res1[1])
anthstd1 = np.nanstd(clim_res1[1])


clim_res2 = np.ones((len(exp2),2))*np.nan
clim_std2 = np.ones((len(exp2),2))*np.nan

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
    respir_calc = dinuptake*(117./16) # too lazy to change var names below
    
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
        clim_std2[e_i,0] = np.nanstd(respir_avg2[((dt>timest3)&(dt<timeen3))])
        clim_std2[e_i,1] = np.nanstd(respir_avg2[((dt>timest4)&(dt<timeen4))])

plt2 = np.nanmean(clim_res2,axis=1)[2:]

anthmean2 = np.nanmean(clim_res2[1])
anthstd2 = np.nanstd(clim_res2[1])

avg98 = clim_res2[2:,0]
avg99 = clim_res2[2:,1]
avg16 = clim_res1[2:,0]
avg17 = clim_res1[2:,1]

std98 = clim_std2[2:,0]
std99 = clim_std2[2:,1]
std16 = clim_std1[2:,0]
std17 = clim_std1[2:,1]

anth98 = clim_res2[1,0]
anth99 = clim_res2[1,1]
anth16 = clim_res1[1,0]
anth17 = clim_res1[1,1]

anthstd98 = clim_std2[1,0]
anthstd99 = clim_std2[1,1]
anthstd16 = clim_std1[1,0]
anthstd17 = clim_std1[1,1]

if sce == 'recy':
    ind1 = list(range(len(title_exp1)))
    ind2 = list(range(len(title_exp2)))
if sce == 'nman':
    ind1 = list(range(len(title_exp1)))
    ind2 = list(range(len(title_exp2)))

if region_name == 'grid':
    yb = -3
    yt = 7
if region_name == 'coast':
    yb = 0
    yt = 25
if region_name == 'offshore':
    yb = -5
    yt = 20

# plot 98
fig,ax = plt.subplots(1,1,figsize=[figw2,figh])

# plot scenarios as bars
ax.bar(ind2,avg98,yerr=std98,color='white',edgecolor='k',capsize=4)
ax.set_title(regtitle+' 1998',fontsize=axfont)

# plot anth as line
ax.plot(ind2,np.ones((len(ind2)))*anth98)
ax.fill_between(ind2,np.ones((len(ind2)))*(anth98+anthstd98),np.ones((len(ind2)))*(anth98-anthstd98),alpha=0.3)

ax.set_ylabel('Uptake mmol C m$^{-3}$ d$^{-1}$',fontsize=axfont)
ax.tick_params(axis='both',which='major',labelsize=axfont)
ax.set_xticks(range(len(title_exp2)))
ax.set_xticklabels(title_exp2,fontsize=axfont,rotation=90)

ax.set_ylim(bottom=yb,top=yt)

savename = 'massb_avg_'+varn+dep+'_uptakeC_'+region_name+'_'+timep+'_'+sce+'98.png'
fig.savefig(savepath+savename,bbox_inches='tight')

# plot 99
fig,ax = plt.subplots(1,1,figsize=[figw2,figh])

# plot scenarios as bars
ax.bar(ind2,avg99,yerr=std99,color='white',edgecolor='k',capsize=4)
ax.set_title(regtitle+' 1999',fontsize=axfont)

# plot anth as line
ax.plot(ind2,np.ones((len(ind2)))*anth99)
ax.fill_between(ind2,np.ones((len(ind2)))*(anth99+anthstd99),np.ones((len(ind2)))*(anth99-anthstd99),alpha=0.3)

ax.set_ylabel('Uptake mmol C m$^{-3}$ d$^{-1}$',fontsize=axfont)
ax.tick_params(axis='both',which='major',labelsize=axfont)
ax.set_xticks(range(len(title_exp2)))
ax.set_xticklabels(title_exp2,fontsize=axfont,rotation=90)

ax.set_ylim(bottom=yb,top=yt)

savename = 'massb_avg_'+varn+dep+'_uptakeC_'+region_name+'_'+timep+'_'+sce+'99.png'
fig.savefig(savepath+savename,bbox_inches='tight')


# plot 16
fig,ax = plt.subplots(1,1,figsize=[figw1,figh])

# plot scenarios as bars
ax.bar(ind1,avg16,yerr=std16,color='gray',capsize=4)
ax.plot(ind1,np.ones((len(ind1)))*anth16)
ax.fill_between(ind1,np.ones((len(ind1)))*(anth16+anthstd16),np.ones((len(ind1)))*(anth16-anthstd16),alpha=0.3)

ax.set_title(regtitle+' 2016',fontsize=axfont)

# plot anth as line
ax.plot(ind1,np.ones((len(ind1)))*anth16)
ax.fill_between(ind1,np.ones((len(ind1)))*(anth16+anthstd16),np.ones((len(ind1)))*(anth16-anthstd16),alpha=0.3)

ax.set_ylabel('Uptake mmol C m$^{-3}$ d$^{-1}$',fontsize=axfont)
ax.tick_params(axis='both',which='major',labelsize=axfont)
ax.set_xticks(range(len(title_exp1)))
ax.set_xticklabels(title_exp1,fontsize=axfont,rotation=90)

ax.set_ylim(bottom=yb,top=yt)

savename = 'massb_avg_'+varn+dep+'_uptakeC_'+region_name+'_'+timep+'_'+sce+'16.png'
fig.savefig(savepath+savename,bbox_inches='tight')

# plot 17
fig,ax = plt.subplots(1,1,figsize=[figw1,figh])

# plot scenarios as bars
ax.bar(ind1,avg17,yerr=std17,color='gray',capsize=4)
ax.plot(ind1,np.ones((len(ind1)))*anth17)
ax.fill_between(ind1,np.ones((len(ind1)))*(anth17+anthstd17),np.ones((len(ind1)))*(anth17-anthstd17),alpha=0.3)

ax.set_title(regtitle+' 2017',fontsize=axfont)

# plot anth as line
ax.plot(ind1,np.ones((len(ind1)))*anth17)
ax.fill_between(ind1,np.ones((len(ind1)))*(anth17+anthstd17),np.ones((len(ind1)))*(anth17-anthstd17),alpha=0.3)

ax.set_ylabel('Uptake mmol C m$^{-3}$ d$^{-1}$',fontsize=axfont)
ax.tick_params(axis='both',which='major',labelsize=axfont)
ax.set_xticks(range(len(title_exp1)))
ax.set_xticklabels(title_exp1,fontsize=axfont,rotation=90)

ax.set_ylim(bottom=yb,top=yt)

savename = 'massb_avg_'+varn+dep+'_uptakeC_'+region_name+'_'+timep+'_'+sce+'17.png'
fig.savefig(savepath+savename,bbox_inches='tight')
