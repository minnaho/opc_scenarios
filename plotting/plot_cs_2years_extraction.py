import sys
sys.path.append('/data/project3/minnaho/global/')
import ROMS_depths as rd
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import glob
import cmocean

#plt.ion()
var_nc = 'O2'
loc = 'bw'

grid_path = '/data/project6/ROMS/L2SCB_AP/roms_grd.nc'
grid_nc = Dataset(grid_path,'r')
pm = np.squeeze(grid_nc['pm'])
pn = np.squeeze(grid_nc['pn'])

sizex = 1E-3/pm
sizey = 1E-3/pn

#catalina island
if loc == 'cat':
    ymin = 450
    ymax = 460

# san nicolas
elif loc == 'sn':
    ymin = 750
    ymax = 760

# bightwide
elif loc == 'bw':
    ymin = 10
    ymax = 1411-10 # remove top and bottom 20 grid points

xplt = np.nanmean(sizex[ymin:ymax,:],axis=0)

cmap1 = cmocean.cm.balance_r

savepath = './figs/2years/'

##################
# cross section
##################
exp = ['cntrl_initap_realistic',
       'fulll_2012_2017',
       'PNDN_only_realistic',
       'pndn50_fixriver',
       'pndn90_fixriver',
       'FNDN_only_realistic',
       'fndn50_fixriver',
       'fndn90_fixriver']
title_exp = ['CTRL',
             'ANTH',
             '50% N Red.',
             '50% N Red.\n50% Recy.',
             '50% N Red.\n90% Recy.',
             '85% N Red.',
             '85% N Red.\n50% Recy.',
             '85% N Red.\n90% Recy.']

# 2 year average
# average y: ymin-ymax, x all points
avg2year_path = '/data/project6/minnaho/opc_scenarios/avg2years/avg_from_extractions/'

axfont=16

################
# plot
# anth vs cntrl
################
if var_nc == 'O2':
    v_mina = -6
    v_maxa = 3
    v_mins = -4
    v_maxs = 4
    depmin = 2
    depmax = 200
if var_nc == 'pH':
    v_mina = -.02
    v_maxa = .02
    v_mins = -.01
    v_maxs = .01
    depmin = 2
    depmax = 200
if var_nc == 'biomass':
    v_mina = -1
    v_maxa = 2.5
    v_mins = -1
    v_maxs = 1
    depmin = 2
    depmax = 100

fig,ax = plt.subplots(1,2,figsize=[12,6])
cntrl_nc = Dataset(avg2year_path+'avg2years_0_200_'+var_nc+'_'+exp[0]+'.nc','r')
fulll_nc = Dataset(avg2year_path+'avg2years_0_200_'+var_nc+'_'+exp[1]+'.nc','r')

depth_arr = np.arange(2,202,2)

data_cntrl = np.squeeze(cntrl_nc['var'])
data_cntrl[data_cntrl>1E10] = np.nan
slice_cntrl = np.nanmean(data_cntrl[:,ymin:ymax,:],axis=1)

data_fulll = np.squeeze(fulll_nc['var'])
data_fulll[data_fulll>1E10] = np.nan
slice_fulll = np.nanmean(data_fulll[:,ymin:ymax,:],axis=1)

pplot0 = ax.flat[0].pcolor((np.arange(data_cntrl.shape[2])-575)*-1*xplt,depth_arr,slice_fulll,cmap=cmocean.cm.ice)
ax.flat[0].set_ylim(top=depmax,bottom=depmin)
ax.flat[0].set_xlim([0,90])
ax.flat[0].invert_xaxis()
ax.flat[0].invert_yaxis()
ax.flat[0].tick_params(axis='both',which='major',labelsize=axfont)
ax.flat[0].set_ylabel('Depth (m)',fontsize=axfont)
ax.flat[0].set_title('ANTH',fontsize=axfont)
p0 = ax.flat[0].get_position().get_points().flatten()
p1 = ax.flat[0].get_position().get_points().flatten()
cb_ax0 = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])
cb0 = fig.colorbar(pplot0,cax=cb_ax0,orientation='vertical')
cb0.ax.tick_params(axis='both',which='major',labelsize=axfont)
cb0.set_label(var_nc+' mmol m$^{-3}$',fontsize=axfont)

#pplot1 = ax.flat[1].pcolor((np.arange(data_fulll.shape[2])-575)*-1*xplt,depth_arr,slice_fulll)
pplot1 = ax.flat[1].pcolor((np.arange(data_fulll.shape[2])-575)*-1*xplt,depth_arr,slice_fulll-slice_cntrl,cmap=cmap1,norm=mcolors.TwoSlopeNorm(vmin=v_mina,vcenter=0,vmax=v_maxa))
#pplot1 = ax.flat[1].pcolor((np.arange(data_fulll.shape[2])-575)*-1*xplt,depth_arr,slice_fulll-slice_cntrl,cmap=cmap1)
ax.flat[1].set_ylim(top=depmax,bottom=depmin)
ax.flat[1].set_xlim([0,90])
ax.flat[1].invert_xaxis()
ax.flat[1].invert_yaxis()
ax.flat[1].tick_params(axis='both',which='major',labelsize=axfont)
ax.flat[1].set_title('ANTH-CTRL',fontsize=axfont)
fig.supxlabel('Distance from coast (km)',fontsize=axfont)
fig.suptitle('2016-2017',fontsize=axfont)

p0 = ax.flat[1].get_position().get_points().flatten()
p1 = ax.flat[1].get_position().get_points().flatten()
cb_ax1 = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])
cb1 = fig.colorbar(pplot1,cax=cb_ax1,orientation='vertical')
cb1.ax.tick_params(axis='both',which='major',labelsize=axfont)
cb1.ax.set_yscale('linear')
cb1.set_label('$\Delta$ '+var_nc+' mmol m$^{-3}$',fontsize=axfont)

ax.flat[1].tick_params(labelleft=False)

plt.savefig(savepath+'cs_'+var_nc+'_anth-cntrl_'+loc+'.png',bbox_inches='tight')

fig,ax = plt.subplots(2,3,sharey=True,sharex=True,figsize=[18,10])
for e_i in range(2,len(exp)):
    avg2year_nc = Dataset(avg2year_path+'avg2years_0_200_'+var_nc+'_'+exp[e_i]+'.nc','r')
    # calculate mean over y=ymin-ymax
    datanc = np.squeeze(avg2year_nc['var'])
    datanc[datanc>1E10] = np.nan
    ncslice = np.nanmean(datanc[:,ymin:ymax,:],axis=1)
    pplot1 = ax.flat[e_i-2].pcolor((np.arange(datanc.shape[2])-575)*-1*xplt,depth_arr,ncslice-slice_fulll,cmap=cmap1,norm=mcolors.TwoSlopeNorm(vmin=v_mins,vcenter=0,vmax=v_maxs))
    #pplot1 = ax.flat[e_i-2].pcolor((np.arange(datanc.shape[2])-575)*-1*xplt,depth_arr,ncslice-slice_fulll,cmap=cmap1)
    ax.flat[e_i-2].set_ylim(top=depmax,bottom=depmin)
    if loc == 'cat':
        ax.flat[e_i-2].set_xlim([0,90])
    if loc == 'sn' or loc == 'bw':
        ax.flat[e_i-2].set_xlim([0,120])
    ax.flat[e_i-2].invert_xaxis()
    ax.flat[e_i-2].invert_yaxis()
    ax.flat[e_i-2].tick_params(axis='both',which='major',labelsize=axfont)
    ax.flat[e_i-2].set_title(title_exp[e_i],fontsize=axfont)

ax.flat[0].set_ylabel('Depth (m)',fontsize=axfont)
ax.flat[3].set_ylabel('Depth (m)',fontsize=axfont)
ax.flat[4].set_xlabel('Distance from coast (km)',fontsize=axfont)
fig.suptitle('Scenario - ANTH',fontsize=axfont)
p0 = ax.flat[2].get_position().get_points().flatten()
p1 = ax.flat[5].get_position().get_points().flatten()
cb_ax1 = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])
cb1 = fig.colorbar(pplot1,cax=cb_ax1,orientation='vertical')
cb1.ax.set_yscale('linear')
cb1.ax.tick_params(axis='both',which='major',labelsize=axfont)
if var_nc != 'pH':
    cb1.set_label('$\Delta$ '+var_nc+' mmol m$^{-3}$',fontsize=axfont)
else:
    cb1.set_label('$\Delta$ '+var_nc,fontsize=axfont)
plt.subplots_adjust(hspace=0.2)
plt.savefig(savepath+'cs_'+var_nc+'_sce-anth_'+loc+'.png',bbox_inches='tight')
plt.close()

# plot profiles
# take average over cross section

lstyle = ['-','--',':','-','--',':']
profcols = ['orange','orange','orange','gray','gray','gray']

fig,ax = plt.subplots(1,2,figsize=[12,6],sharey=True)
for e_i in range(2,len(exp)):
    avg2year_nc = Dataset(avg2year_path+'avg2years_0_200_'+var_nc+'_'+exp[e_i]+'.nc','r')
    # calculate mean over y=ymin-ymax
    datanc = np.squeeze(avg2year_nc['var'])
    datanc[datanc>1E10] = np.nan
    #ncslice = np.nanmean(datanc[:,ymin:ymax,:],axis=1)
    #prof_exp = np.nanmean(ncslice-slice_fulll,axis=1)
    prof_exp = np.nanmean(datanc-data_fulll,axis=(1,2))
    
    if e_i in range(2,5):
        ax.flat[0].plot(prof_exp,depth_arr,linestyle=lstyle[e_i-2],linewidth=3,color=profcols[e_i-2],label=title_exp[e_i])
    else: 
        ax.flat[1].plot(prof_exp,depth_arr,linestyle=lstyle[e_i-2],linewidth=3,color=profcols[e_i-2],label=title_exp[e_i])

if var_nc == 'O2' or var_nc == 'biomass':
    #ax.flat[0].set_xlabel('$\Delta$ '+var_nc+' mmol m$^{-3}$',fontsize=axfont)
    fig.supxlabel('$\Delta$ '+var_nc+' mmol m$^{-3}$',fontsize=axfont)

if var_nc == 'pH':
    #ax.flat[0].set_xlabel('$\Delta$ '+var_nc,fontsize=axfont)
    fig.supxlabel('$\Delta$ '+var_nc,fontsize=axfont)
    ax.flat[0].set_xlim([-0.0015,0.0025])
    ax.flat[1].set_xlim([-0.0015,0.0025])

if var_nc == 'O2':
    if loc == 'cat':
        ax.flat[0].set_xlim([-1.5,3])
        ax.flat[1].set_xlim([-1.5,3])
    if loc == 'sn':
        ax.flat[0].set_xlim([-3.5,2.5])
        ax.flat[1].set_xlim([-3.5,2.5])
    if loc == 'bw':
        ax.flat[0].set_xlim([-0.75,1.5])
        ax.flat[1].set_xlim([-0.75,1.5])

if var_nc == 'biomass':
    ax.flat[0].set_xlim([-0.5,0.6])
    ax.flat[0].set_ylim([0,80])
    ax.flat[1].set_xlim([-0.5,0.6])
    ax.flat[1].set_ylim([0,80])

ax.flat[0].invert_yaxis()
#ax.flat[1].invert_yaxis()

ax.flat[0].set_ylabel('Depth (m)',fontsize=axfont)
#ax.set_title('Scenario - ANTH',fontsize=axfont)
ax.flat[0].tick_params(axis='both',which='major',labelsize=axfont)
ax.flat[1].tick_params(axis='both',which='major',labelsize=axfont)
if var_nc == 'pH':
    ax.flat[0].legend(loc='best',fontsize=axfont-2)
    #ax.flat[0].set_xticklabels(ax.get_xticklabels(),rotation=60,ha='right')
    ax.flat[1].legend(loc='best',fontsize=axfont-2)
    #ax.flat[1].set_xticklabels(ax.get_xticklabels(),rotation=60,ha='right')

if var_nc == 'O2':
    ax.flat[0].legend(loc='lower right',fontsize=axfont-2)
    ax.flat[1].legend(loc='upper right',fontsize=axfont-2)

fig.savefig(savepath+'cs_profile_'+var_nc+'_'+loc+'_split_newmean.png',bbox_inches='tight')
plt.close()

