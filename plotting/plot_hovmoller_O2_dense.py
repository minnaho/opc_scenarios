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
var = 'rho'
var_nc = 'rho'
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

figpath = './figs/hovmollers/'


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

difexp = np.load(savenpypath+'dif_hov_'+var+'_'+exp[-1]+'_'+region_name+'.npy')
absexp = np.load(savenpypath+'abs_hov_'+var+'_'+exp[-1]+'_'+region_name+'.npy')
abscnt = np.load(savenpypath+'abs_hov_'+var+'_'+substr+'_'+region_name+'.npy')

time_units = 'days since 1997-11-01'

dtplt = num2date(range(difexp.shape[1]),time_units,only_use_cftime_datetimes=False,only_use_python_datetimes=True)

rho0 = 1024.7


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
    cblabel = var+' mmol change'
if var == 'omega':
    v_max = 0.2
    v_min = -0.2
    cblabel = var+' change'
if var == 'pH':
    v_max = 0.05
    v_min = -0.05
    cblabel = var+' change'

# multiply by -1 for now because switched the subtraction in 
# extract hovmoller script
plt_trans = np.transpose(difexp,axes=[0,2,1])*-1
plt_absexp = np.transpose(absexp,axes=[0,2,1])
plt_abscnt = np.transpose(abscnt,axes=[0,2,1])

# make mask to contour over
val = 1.2
contour = np.array((np.where((plt_absexp[0]<val)&(plt_abscnt[0]>val))[0],np.where((plt_absexp[0]<val)&(plt_abscnt[0]>val))[1]))
mask = np.ones((plt_trans.shape[1],plt_trans.shape[2]))
mask[contour[0],contour[1]] = 0

depplt = range(2,202,2)

if len(exp) > 1:
    figw = 14
    figh = 10
    fig,ax = plt.subplots(len(exp),1,figsize=[figw,figh])
    for n_i in range(len(exp)):
        p_plt = ax[n_i].pcolormesh(dtplt,depplt,plt_trans[n_i],cmap=c_map,vmin=v_min,vmax=v_max)
        ax[n_i].set_title(title_exp[n_i]+' - CTRl',fontsize=axisfont)
        #ax[n_i].plot(dtplt,np.ones(roms_plt.shape[0])*0,color='k',linestyle='--')
        #ax[n_i].set_ylim([v_min,v_max])
        #if perc == False:
        #    tick_spacingy = 100
        #else:
        #    tick_spacingy = 30
        #ax[n_i].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacingy))
    
        ax[n_i].invert_yaxis()
        ax[n_i].set_ylabel(ylabel,fontsize=axisfont)
        ax[n_i].tick_params(axis='both',which='major',labelsize=axisfont)
        ax[n_i].yaxis.set_ticks_position('both')
        ax[n_i].xaxis.set_ticks_position('both')
        if n_i < len(exp)-1:
            ax[n_i].xaxis.set_ticklabels([])
        p0 = ax.flat[0].get_position().get_points().flatten()
        p1 = ax.flat[-1].get_position().get_points().flatten()
        cb_ax = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])

        cb = fig.colorbar(p_plt,cax=cb_ax,orientation='vertical')
        cb.set_label(cblabel,fontsize=axisfont)
        cb.ax.tick_params(axis='both',which='major',labelsize=axisfont)
else:
    figw = 14
    figh = 5
    fig,ax = plt.subplots(len(exp),1,figsize=[figw,figh])
    # hovmoller
    p_plt = ax.pcolormesh(dtplt,depplt,plt_trans[0],cmap=c_map,vmin=v_min,vmax=v_max)
    # contour
    c_plt = ax.contour(dtplt,depplt,mask,[0])
    ax.set_title(title_exp[0]+' - CTRL',fontsize=axisfont)

    ax.invert_yaxis()
    ax.set_ylabel(ylabel,fontsize=axisfont)
    ax.tick_params(axis='both',which='major',labelsize=axisfont)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    p0 = ax.get_position().get_points().flatten()
    p1 = ax.get_position().get_points().flatten()
    cb_ax = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])

    cb = fig.colorbar(p_plt,cax=cb_ax,orientation='vertical')
    cb.set_label(cblabel,fontsize=axisfont)
    cb.ax.tick_params(axis='both',which='major',labelsize=axisfont)


savename = 'hov_dif_'+var+'_'+region_name+'_'+exp[-1]+'_contour.png'
fig.savefig(figpath+savename,bbox_inches='tight')
