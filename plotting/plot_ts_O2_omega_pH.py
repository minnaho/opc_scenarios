# plot timeseries of integrated and
# sliced varibles
import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import l2grid
import numpy as np
from netCDF4 import Dataset,num2date
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import glob as glob


plt.ion()

region_name = 'grid'
var = 'O2'
#var = 'omega_pH'
var_nc = 'var'
#var_nc = 'pH'
dpst = '0'
dpen = '200'

ncpath = '/data/project6/minnaho/opc_scenarios/ext_depth_200/'
filest = 'ext_'+dpst+'_'+dpen+'_'+var+'_'

fileen = '.nc'

# experiment to subtract from
substr = 'cntrl'

figpath = './figs/ts/'

ylabel = var+' mmol loss'

# choose years
start_year = 1997
end_year = 1999

# choose months between 1 and 12
start_month = 11
end_month = 11 

# scenario names
#exp = ['l1617','PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
#title_exp = ['Loads 16-17','PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']
#exp = ['PNDN_only','pndn50','pndn90']
#title_exp = ['PNDN only','PNDN 50','PNDN 90']
#exp = ['FNDN_only','fndn50','fndn90']
#title_exp = ['FNDN only','FNDN 50','FNDN 90']
exp = ['l1617']
title_exp = ['Loads 16-17']

time_units = 'days since 1997-11-01'

# last day of simulation to create time dimension of plot
timenc = np.array(Dataset(ncpath+filest+'Y'+str(end_year)+'M'+str(end_month)+'D27_l1617.nc','r').variables['time'])
dtplt = num2date(range(timenc[0].astype(int)+1),time_units,only_use_cftime_datetimes=False,only_use_python_datetimes=True)



mask_nc = l2grid.mask_nc
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

pltnc = np.ones((len(exp),len(dtplt)))*np.nan
t_i = 0
for y_i in range(start_year,end_year+1):
    # if we are on the first year, starts at s_m
    if y_i == start_year:
        s_m = start_month
    else:
        s_m = 1
    # if we are on the last year, end at e_m
    if y_i == end_year:
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
        if y_i == 1999 and m_i == 11:
            ndays = 27 # last day of the simulation
        for d_i in list(range(1,ndays+1)):
            dtstr = 'Y'+str(y_i)+'M'+'%02d'%m_i+'D'+'%02d'%d_i
            print(exp,dtstr)
            for n_i in range(len(exp)):
                expnc = np.squeeze(Dataset(ncpath+filest+dtstr+'_'+exp[n_i]+fileen,'r').variables[var_nc])*mask_mult
                subnc = np.squeeze(Dataset(ncpath+filest+dtstr+'_'+substr+fileen,'r').variables[var_nc])*mask_mult
                diffnc = expnc-subnc
                # sum over negative values then make the value positive
                pltnc[n_i][t_i] = np.nansum(diffnc[diffnc<0])*-1
            t_i += 1

if var == 'O2':
    np.save('ts_sum_neg_'+var+'_'+exp[-1]+'.npy',pltnc)
else:
    np.save('ts_sum_neg_'+var_nc+'_'+exp[-1]+'.npy',pltnc)
    

# plotting
axisfont = 16
figw = 14
figh = 10

fig,ax = plt.subplots(len(exp),1,figsize=[figw,figh])

for n_i in range(len(exp)):

    ax[n_i].plot(dtplt,pltnc[n_i],color='k')
    ax[n_i].set_title(title_exp[n_i],fontsize=axisfont)
    #ax[n_i].plot(dtplt,np.ones(roms_plt.shape[0])*0,color='k',linestyle='--')
    #ax[n_i].set_ylim([v_min,v_max])
    #if perc == False:
    #    tick_spacingy = 100
    #else:
    #    tick_spacingy = 30
    #ax[n_i].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacingy))

    ax[n_i].set_ylabel(ylabel,fontsize=axisfont)
    ax[n_i].tick_params(axis='both',which='major',labelsize=axisfont)
    ax[n_i].yaxis.set_ticks_position('both')
    ax[n_i].xaxis.set_ticks_position('both')
    if n_i < len(exp)-1:
        ax[n_i].xaxis.set_ticklabels([])


savename = 'ts_sum_neg_'+var+'_'+dp+'_'+region_name+'_'+exp[n_i]+'.png'
fig.savefig(figpath+savename,bbox_inches='tight')
