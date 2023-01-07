import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import l2grid
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import calendar

plt.ion()

outpath = '/data/project6/minnaho/opc_scenarios/pdf_npy/'
savepath = './figs/pdf/'

region_name = 'coast'

varstr1 = 'dens'
varstr2 = 'salt'

start_year1 = 1998
end_year1 = 1998

start_year2 = 2016
end_year2 = 2016

start_month = 3
end_month = 3

savename1 = 'pdf_'+varstr1+'_'+str(start_year1)+'_'+str(start_year2)+'M'+'%02d'%start_month
savename2 = 'pdf_'+varstr2+'_'+str(start_year1)+'_'+str(start_year2)+'M'+'%02d'%start_month

exp1 = ['cntrl','loads1617','PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
exp2 = ['cntrl_2012_2017','fulll_2012_2017','PNDN_only_realistic','pndn50_realistic','pndn90_realistic','FNDN_only_realistic','fndn50_realistic','fndn90_realistic']
title_exp = ['Ocean Only','Current Day','50% N Red.','50% N Red. 50% Recy.','50% N Red. 90% Recy.','85% N Red.','85% N Red. 50% Recy.','85% N Red. 90% Recy.']

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
mask_nc = l2grid.mask_nc # full grid

mask_ssd[mask_ssd==0] = np.nan
mask_nsd[mask_nsd==0] = np.nan
mask_oc[mask_oc==0] = np.nan
mask_sp[mask_sp==0] = np.nan
mask_sm[mask_sm==0] = np.nan
mask_v[mask_v==0] = np.nan
mask_sb[mask_sb==0] = np.nan
mask_nc[mask_nc==0] = np.nan

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
if region_name == 'grid': # full L2 grid
    mask_mult = mask_nc
    regtitle = 'SCB'
if region_name == 'coast':
    mask_mult = mask_cst
    regtitle = '15 km coast'

months_w_31_days = [1,3,5,7,8,10,12]
leap_years = [1992,1996,2000,2004,2008,2012,2016,2020]

nbinshape = 501

for y_i in range(start_year1,end_year1+1):
    # if we are on the first year, starts at s_m
    if y_i == start_year1:
        s_m = start_month
    else:
        s_m = 1
    # if we are on the last year, end at e_m
    if y_i == end_year1:
        e_m = end_month+1
    else:
        e_m = 13
    for m_i in range(s_m,e_m):
        # loop through each file type
        if m_i in months_w_31_days:
            ndays = 31
        if m_i not in months_w_31_days:
            ndays = 30
            if m_i == 2 and y_i in leap_years:
                ndays = 29
            if m_i == 2 and y_i not in leap_years:
                ndays = 28
        n_p_salt1 = np.ones((len(exp1),ndays,nbinshape))*np.nan
        n_p_rho1 = np.ones((len(exp1),ndays,nbinshape))*np.nan
        flt_salt1 = np.ones((len(exp1),ndays,nbinshape))*np.nan
        flt_rho1 = np.ones((len(exp1),ndays,nbinshape))*np.nan
        bin_p_salt1 = np.ones((len(exp1),ndays,nbinshape))*np.nan
        bin_p_rho1 = np.ones((len(exp1),ndays,nbinshape))*np.nan
        for d_i in list(range(1,ndays+1)):
            dt = 'Y'+str(y_i)+'M'+'%02d'%m_i+'D'+'%02d'%d_i
            print(dt)
            for e_i in range(len(exp1)):
                n_p_salt1[e_i,d_i-1,:] = np.load(outpath+'n_p_count_salt_'+dt+'_100m_'+exp1[e_i]+'.npy')
                n_p_rho1[e_i,d_i-1,:]  = np.load(outpath+'n_p_count_rho_'+dt+'_100m_'+exp1[e_i]+'.npy')
                flt_salt1[e_i,d_i-1,:] = np.load(outpath+'flt_nonan_salt_'+dt+'_100m_'+exp1[e_i]+'.npy')
                flt_rho1[e_i,d_i-1,:]  = np.load(outpath+'flt_nonan_rho_'+dt+'_100m_'+exp1[e_i]+'.npy')
                bin_p_salt1[e_i,d_i-1,:] = np.load(outpath+'bin_p_salt_'+dt+'_100m_'+exp1[e_i]+'.npy')
                bin_p_rho1[e_i,d_i-1,:]  = np.load(outpath+'bin_p_rho_'+dt+'_100m_'+exp1[e_i]+'.npy')

for y_i in range(start_year2,end_year2+1):
    # if we are on the first year, starts at s_m
    if y_i == start_year2:
        s_m = start_month
    else:
        s_m = 1
    # if we are on the last year, end at e_m
    if y_i == end_year2:
        e_m = end_month+1
    else:
        e_m = 13
    for m_i in range(s_m,e_m):
        # loop through each file type
        if m_i in months_w_31_days:
            ndays = 31
        if m_i not in months_w_31_days:
            ndays = 30
            if m_i == 2 and y_i in leap_years:
                ndays = 29
            if m_i == 2 and y_i not in leap_years:
                ndays = 28
        n_p_salt2 = np.ones((len(exp2),ndays,nbinshape))*np.nan
        n_p_rho2 = np.ones((len(exp2),ndays,nbinshape))*np.nan
        flt_salt2 = np.ones((len(exp2),ndays,nbinshape))*np.nan
        flt_rho2 = np.ones((len(exp2),ndays,nbinshape))*np.nan
        bin_p_salt2 = np.ones((len(exp2),ndays,nbinshape))*np.nan
        bin_p_rho2 = np.ones((len(exp2),ndays,nbinshape))*np.nan
        for d_i in list(range(1,ndays+1)):
            dt = 'Y'+str(y_i)+'M'+'%02d'%m_i+'D'+'%02d'%d_i
            print(dt)
            for e_i in range(len(exp2)):
                n_p_salt2[e_i,d_i-1,:] = np.load(outpath+'n_p_count_salt_'+dt+'_100m_'+exp2[e_i]+'.npy')
                n_p_rho2[e_i,d_i-1,:]  = np.load(outpath+'n_p_count_rho_'+dt+'_100m_'+exp2[e_i]+'.npy')
                flt_salt2[e_i,d_i-1,:] = np.load(outpath+'flt_nonan_salt_'+dt+'_100m_'+exp2[e_i]+'.npy')
                flt_rho2[e_i,d_i-1,:]  = np.load(outpath+'flt_nonan_rho_'+dt+'_100m_'+exp2[e_i]+'.npy')
                bin_p_salt2[e_i,d_i-1,:] = np.load(outpath+'bin_p_salt_'+dt+'_100m_'+exp2[e_i]+'.npy')
                bin_p_rho2[e_i,d_i-1,:]  = np.load(outpath+'bin_p_rho_'+dt+'_100m_'+exp2[e_i]+'.npy')

n_p_concat_salt = np.nansum(n_p_salt1+n_p_salt2,axis=1)
n_p_concat_rho = np.nansum(n_p_rho1+n_p_rho2,axis=1)

flt_concat_salt = np.nansum(flt_salt1+flt_salt2,axis=1)
flt_concat_rho  = np.nansum(flt_rho1+flt_rho2,axis=1)

# bins are the same between the two yeras and at all time steps
bin_p_salt = bin_p_salt1[0]
bin_p_rho = bin_p_rho1[0]


figw = 12
figh = 8
axisfont = 16
cplt = ['green','blue','orange','orange','orange','gray','gray','gray']
lsty = ['-','-','-','--',':','-','--',':']

fig,ax = plt.subplots(1,1,figsize=[figw,figh])
for e_i in range(len(title_exp)):
    ax.plot(bin_p_rho[e_i],n_p_concat_rho[e_i]/flt_concat_rho[e_i],color=cplt[e_i],linestyle=lsty[e_i],linewidth=1.5,label=title_exp[e_i])

ax.legend(fontsize=axisfont)
ax.set_xlabel('Density kg/m3',fontsize=axisfont)
ax.set_ylabel('PDF',fontsize=axisfont)
ax.set_title(regtitle+' Density PDF N = 2 years '+calendar.month_abbr[start_month]+' '+str(start_year1)+' '+str(start_year2),fontsize=axisfont)

ax.tick_params(axis='both',which='major',labelsize=axisfont)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
#ax.set_xlim([90,300])
ax.set_yscale('log')
ax.set_ylim([1E-6,1E-1])
ax.set_xlim([1022,1024.5])

plt.savefig(savepath+savename1+'.png',bbox_inches='tight')

figw = 12
figh = 8
axisfont = 16
cplt = ['green','blue','orange','orange','orange','gray','gray','gray']
lsty = ['-','-','-','--',':','-','--',':']

fig,ax = plt.subplots(1,1,figsize=[figw,figh])
for e_i in range(len(title_exp)):
    ax.plot(bin_p_salt[e_i],n_p_concat_salt[e_i]/flt_concat_salt[e_i],color=cplt[e_i],linestyle=lsty[e_i],linewidth=1.5,label=title_exp[e_i])

ax.legend(fontsize=axisfont)
ax.set_xlabel('Salinity PSU',fontsize=axisfont)
ax.set_ylabel('PDF',fontsize=axisfont)
ax.set_title(regtitle+' Salinity PDF N = 2 years '+calendar.month_abbr[start_month]+' '+str(start_year1)+' '+str(start_year2),fontsize=axisfont)

ax.tick_params(axis='both',which='major',labelsize=axisfont)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
#ax.set_xlim([90,300])
ax.set_yscale('log')
ax.set_ylim([1E-6,1E-1])
ax.set_xlim([30,33])

plt.savefig(savepath+savename2+'.png',bbox_inches='tight')
