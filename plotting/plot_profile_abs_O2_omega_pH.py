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
import cmocean


plt.ion()

region_name = 'grid'
var = 'omega'
var_nc = 'var'
#var = 'omega_pH'
#var_nc = 'omega'

dpst = '0'
dpen = '200'

c_map = 'bwr'

ncpath = '/data/project6/minnaho/opc_scenarios/ext_depth_200/'
filest = 'ext_'+dpst+'_'+dpen+'_'+var+'_'

fileen = '.nc'


# experiment to subtract from
substr = 'cntrl'

figpath = './figs/profiles/'


ylabel = 'Depth (m)'

# choose years
start_year = 1997
end_year = 1999

# choose months between 1 and 12
start_month = 11
end_month = 11 
# scenario names
#exp = ['l1617','PNDN_only','pndn50','pndn90','FNDN_only','fndn50','fndn90']
#title_exp = ['Loads 16-17','PNDN only','PNDN 50','PNDN 90','FNDN only','FNDN 50','FNDN 90']
#exp = ['l1617','PNDN_only','pndn50','pndn90']
#title_exp = ['Loads 16-17','PNDN only','PNDN 50','PNDN 90']
#exp = ['FNDN_only','fndn50','fndn90']
#title_exp = ['FNDN only','FNDN 50','FNDN 90']
#exp = ['PNDN_only','pndn50','pndn90']
#title_exp = ['PNDN only','PNDN 50','PNDN 90']
exp = ['l1617']
title_exp = ['Loads 16-17']

time_units = 'days since 1997-11-01'

savenpypath = '/data/project6/minnaho/opc_scenarios/hovmollers/'

#difexp = np.load(savenpypath+'dif_hov_'+var+'_'+exp[-1]+'_'+region_name+'.npy')
absexp = np.load(savenpypath+'abs_hov_'+var+'_'+exp[-1]+'_'+region_name+'.npy')
abscnt = np.load(savenpypath+'abs_hov_'+var+'_'+substr+'_'+region_name+'.npy')

# plotting
axisfont = 16


#if perc == False:
#    v_max = 220
#    v_min = -130
#else:
#    v_max = 170
#    v_min = -45

if var == 'O2':
    v_max = 30
    v_min = -30
    cblabel = var+' mmol'
if var == 'omega':
    v_max = 0.2
    v_min = -0.2
    cblabel = var
if var == 'pH':
    v_max = 0.05
    v_min = -0.05
    cblabel = var

 
plt_trans = np.transpose(absexp,axes=[0,2,1])
plt_prof = np.nanmean(plt_trans,axis=2)
plt_prof_05 = np.nanpercentile(plt_trans,5,axis=2)
plt_prof_95 = np.nanpercentile(plt_trans,95,axis=2)

plt_trans_cnt = np.transpose(abscnt,axes=[0,2,1])
plt_prof_cnt = np.nanmean(plt_trans_cnt,axis=2)
plt_prof_05_cnt = np.nanpercentile(plt_trans_cnt,5,axis=2)
plt_prof_95_cnt = np.nanpercentile(plt_trans_cnt,95,axis=2)

depplt = range(2,202,2)

linee = 'cornflowerblue'
linec = 'orange'

figw = 7
figh = 12

fig,ax = plt.subplots(len(exp),1,figsize=[figw,figh])
ax.plot(plt_prof[0],depplt,color=linee)
ax.plot(plt_prof_05[0],depplt,color=linee,linestyle='--')
ax.plot(plt_prof_95[0],depplt,color=linee,linestyle='--')
ax.plot(plt_prof_cnt[0],depplt,color=linec)
ax.plot(plt_prof_05_cnt[0],depplt,color=linec,linestyle='--')
ax.plot(plt_prof_95_cnt[0],depplt,color=linec,linestyle='--')
ax.fill_betweenx(depplt,plt_prof_05[0],plt_prof_95[0],color=linec,alpha=.3)
ax.fill_betweenx(depplt,plt_prof_05_cnt[0],plt_prof_95_cnt[0],color=linec,alpha=.3)
ax.set_title(title_exp[0]+' - CTRL',fontsize=axisfont)
ax.invert_yaxis()
ax.set_ylabel(ylabel,fontsize=axisfont)
ax.tick_params(axis='both',which='major',labelsize=axisfont)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.set_xlabel(cblabel,fontsize=axisfont)


savename = 'prof_abs_'+var+'_'+region_name+'_'+exp[-1]+'.png'
fig.savefig(figpath+savename,bbox_inches='tight')
