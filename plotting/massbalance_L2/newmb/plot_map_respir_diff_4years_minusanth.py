# change to crtopy_update environment before running
# conda activate cartopy_update
# N = 4 years
# compare to ANTH
import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import l2grid as l2grid
import numpy as np
import h5py
from netCDF4 import Dataset
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cpf
import cmocean

#plt.ion()

region_name = 'grid'

outpath = './budget_ww/'
savepath = './figs/'

timep = 'fullts'

sce = 'recy'

c_map = cmocean.cm.balance

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
    # N = 4
    exp_name = ['ANTH',
           'CTRL',
           'PNDN_only',
           'pndn50',
           'pndn90',
           'FNDN_only',
           'fndn50',
           'fndn90'
                                ]
    exp1 = [anth1,
           cntrl1,
           'PNDN_only_realistic',
           'pndn50_realistic',
           'pndn90_realistic',
           'FNDN_only_realistic',
            'fndn50_realistic',
            'fndn90_realistic']

    exp2 = [anth2,
           cntrl2,
           'PNDN_only',
           'pndn50',
           'pndn90',
           'FNDN_only',
           'fndn50',
           'fndn90'
                ]
    title_exp1 = ['CTRL',
                 '50% N Red.',
                 '50% N Red.\n50% Recy.',
                 '50% N Red.\n90% Recy.',
                 '85% N Red.',
                 '85% N Red.\n50% Recy.',
                 '85% N Red.\n90% Recy.'
                ]
    title_exp2 = ['Land-based',
                 '50% N Red.',
                 '50% N Red.\n50% Recy.',
                 '50% N Red.\n90% Recy.',
                 '85% N Red.',
                 '85% N Red.\n50% Recy.',
                 '85% N Red.\n90% Recy.'
                ]
                            
    #title_exp = [
    #             '50% N Red.',
    #             '50% N Red.\n50% Recy.',
    #             '50% N Red.\n90% Recy.',
    #             '85% N Red.',
    #             '85% N Red.\n50% Recy.',
    #             '85% N Red.\n90% Recy.',
    #             '50% N Red.',
    #             '50% N Red.\n50% Recy.',
    #             '50% N Red.\n90% Recy.',
    #             '85% N Red.',
    #            ]

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


varn = 'O2'
matn = 'MATBGCF'
matc = 'MATVARC'

fst = 'outputs_'+varn+'_'
dep = '_0-200m'


mask_nc = l2grid.mask_nc

s2d = 86400

figw = 10
figh = 8
axfont = 14

# calculate BGC and uptake terms
clim_res1 = np.ones((len(exp1),2,mask_nc.shape[0],mask_nc.shape[1]))*np.nan

for e_i in range(len(exp1)):
    # read in file
    matr = h5py.File(outpath+fst+exp1[e_i]+dep+'/'+matn+'.mat','r')
    data = matr.get(matn)

    # get dates
    datemat = np.squeeze(h5py.File(outpath+fst+exp1[e_i]+dep+'/'+matc+'.mat','r').get(matc)['date'])

    # get variables
    loss = np.squeeze(data['LOSS'])
    graze = np.squeeze(data['GRAZE'])
    remin = np.squeeze(data['REMIN'])
    sedre = np.squeeze(data['SED_REMIN'])
    ammox = np.squeeze(data['AMMOX'])
    nit = np.squeeze(data['NIT'])

    # calculate terms
    respir_calc = loss+graze+remin+sedre+ammox+nit
    
    # subtract from experiment
    if exp1[e_i] == anth1:
        respir_cntrl1 = np.copy(respir_calc)
        respir_cntrl1[respir_cntrl1==0] = np.nan

    else:
        respir_avg1 = (respir_calc-respir_cntrl1)*s2d

        respir_avg1[respir_avg1==0] = np.nan

        # convert matlab time
        dt = pd.to_datetime(datemat-719529, unit='D')

        # calculate avg and std over each year
        clim_res1[e_i,0] = np.nanmean(respir_avg1[((dt>timest1)&(dt<timeen1))],axis=0)
        clim_res1[e_i,1] = np.nanmean(respir_avg1[((dt>timest2)&(dt<timeen2))],axis=0)


clim_res2 = np.ones((len(exp2),2,mask_nc.shape[0],mask_nc.shape[1]))*np.nan

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
    if exp2[e_i] == anth2:
        respir_cntrl2 = np.copy(respir_calc)
        respir_cntrl2[respir_cntrl2==0] = np.nan

    else:
    
        respir_avg2 = (respir_calc-respir_cntrl2)*s2d

        respir_avg2[respir_avg2==0] = np.nan

        # convert matlab time
        dt = pd.to_datetime(datemat-719529, unit='D')

        # calculate avg and std over each year
        clim_res2[e_i,0] = np.nanmean(respir_avg2[((dt>timest3)&(dt<timeen3))],axis=0)
        clim_res2[e_i,1] = np.nanmean(respir_avg2[((dt>timest4)&(dt<timeen4))],axis=0)

clim_res1[clim_res1>1E10] = np.nan
clim_res2[clim_res2>1E10] = np.nan

clim_res1[clim_res1==0] = np.nan
clim_res2[clim_res2==0] = np.nan

# plot difference 
plt1 = clim_res1[1:]
plt2 = clim_res2[1:]

# take N = 4 mean for each exp
avg4year = np.ones((clim_res1[1:].shape[0],mask_nc.shape[0],mask_nc.shape[1]))
for c_i in range(1,clim_res1.shape[0]):
    avg4year[c_i-1] = np.nanmean(np.concatenate((clim_res1[c_i],clim_res2[c_i]),axis=0),axis=0)

latnc = l2grid.lat_nc
lonnc = l2grid.lon_nc

# full grid
if region_name == 'grid':
    lat_min = 32.5
    lat_max = 35
    lon_min = -121
    lon_max = -117

extent = [lon_min,lon_max,lat_min,lat_max]
coast_10m = cpf.NaturalEarthFeature('physical','coastline','10m')

# plot

#cblabel = 'Respiration mmol O m$^{-3}$ d$^{-1}$'
cblabel = 'Respiration mmol O m$^{-3}$ d$^{-1}$'

v_min = -20
v_max = 20

for e_p in range(len(title_exp1)):
    fig,ax = plt.subplots(1,1,figsize=[figw,figh],subplot_kw={'projection':ccrs.PlateCarree()})
    pplot = ax.pcolormesh(lonnc,latnc,avg4year[e_p],transform=ccrs.PlateCarree(),cmap=c_map,vmin=v_min,vmax=v_max)
    ax.add_feature(coast_10m,facecolor='None',edgecolor='k')
    ax.add_feature(cpf.BORDERS,facecolor='None',edgecolor='k')
    ax.set_extent(extent)
    gl = ax.gridlines(draw_labels={'bottom':'x','left':'y'},xlabel_style=dict(size=axfont),ylabel_style=dict(size=axfont))
    ax.set_title(title_exp1[e_p],fontsize=axfont)
    p0 = ax.get_position().get_points().flatten()
    p1 = ax.get_position().get_points().flatten()
    cb_ax = fig.add_axes([p0[2]+.01,p1[1],.01,p0[3]-p1[1]])
    cb = fig.colorbar(pplot,cax=cb_ax,orientation='vertical')
    cb.set_label(cblabel,fontsize=axfont)
    cb.ax.tick_params(axis='both',which='major',labelsize=axfont)
    savename = 'map_4years_minusanth_massb_diff_'+varn+dep+'_'+region_name+'_'+timep+'_'+exp_name[e_p+1]+'.png'
    fig.savefig(savepath+savename,bbox_inches='tight')

'''
for e_p in range(len(title_exp1)):
    for y_p in range(clim_res1.shape[1]):
        fig,ax = plt.subplots(1,1,figsize=[figw,figh],subplot_kw={'projection':ccrs.PlateCarree()})
        pplot = ax.pcolormesh(lonnc,latnc,plt1[e_p,y_p],transform=ccrs.PlateCarree(),cmap=c_map,vmin=v_min,vmax=v_max)
        ax.add_feature(coast_10m,facecolor='None',edgecolor='k')
        ax.add_feature(cpf.BORDERS,facecolor='None',edgecolor='k')
        ax.set_extent(extent)
        gl = ax.gridlines(draw_labels={'bottom':'x','left':'y'},xlabel_style=dict(size=axfont),ylabel_style=dict(size=axfont))
        if y_p == 0:
            yr = str(2016)
        if y_p == 1:
            yr = str(2017)
        ax.set_title(title_exp1[e_p]+' '+yr,fontsize=axfont)
        p0 = ax.get_position().get_points().flatten()
        p1 = ax.get_position().get_points().flatten()
        cb_ax = fig.add_axes([p0[2]+.01,p1[1],.01,p0[3]-p1[1]])
        cb = fig.colorbar(pplot,cax=cb_ax,orientation='vertical')
        cb.set_label(cblabel,fontsize=axfont)
        cb.ax.tick_params(axis='both',which='major',labelsize=axfont)
        savename = 'map_massb_diff_'+varn+dep+'_'+region_name+'_'+timep+'_'+exp1[e_p+1]+yr+'.png'
        fig.savefig(savepath+savename,bbox_inches='tight')

for e_p in range(len(title_exp2)):
    for y_p in range(clim_res2.shape[1]):
        fig,ax = plt.subplots(1,1,figsize=[figw,figh],subplot_kw={'projection':ccrs.PlateCarree()})
        pplot = ax.pcolormesh(lonnc,latnc,plt2[e_p,y_p],transform=ccrs.PlateCarree(),cmap=c_map,vmin=v_min,vmax=v_max)
        ax.add_feature(coast_10m,facecolor='None',edgecolor='k')
        ax.add_feature(cpf.BORDERS,facecolor='None',edgecolor='k')
        ax.set_extent(extent)
        gl = ax.gridlines(draw_labels={'bottom':'x','left':'y'},xlabel_style=dict(size=axfont),ylabel_style=dict(size=axfont))
        if y_p == 0:
            yr = str(1998)
        if y_p == 1:
            yr = str(1999)
        ax.set_title(title_exp2[e_p]+' '+yr,fontsize=axfont)
        p0 = ax.get_position().get_points().flatten()
        p1 = ax.get_position().get_points().flatten()
        cb_ax = fig.add_axes([p0[2]+.01,p1[1],.01,p0[3]-p1[1]])
        cb = fig.colorbar(pplot,cax=cb_ax,orientation='vertical')
        cb.set_label(cblabel,fontsize=axfont)
        cb.ax.tick_params(axis='both',which='major',labelsize=axfont)
        savename = 'map_massb_diff_'+varn+dep+'_'+region_name+'_'+timep+'_'+exp2[e_p+1]+yr+'.png'
        fig.savefig(savepath+savename,bbox_inches='tight')
'''
