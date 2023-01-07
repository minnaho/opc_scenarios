import sys
sys.path.append('/data/project3/minnaho/global/')
import ROMS_depths as rd
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import cmocean

plt.ion()
var_nc = 'biomass'

grid_path = '/data/project6/ROMS/L2SCB_AP/roms_grd.nc'
grid_nc = Dataset(grid_path,'r')
pm = np.squeeze(grid_nc['pm'])
pn = np.squeeze(grid_nc['pn'])

sizex = 1E-3/pm
sizey = 1E-3/pn

xplt = np.nanmean(sizex[450:460,:],axis=0)

cmap1 = cmocean.cm.balance

##################
# cross section
##################
exp = ['cntrl_initap_realistic','fulll_2012_2017','PNDN_only_realistic','pndn50_realistic','pndn90_realistic','FNDN_only_realistic','fndn50_realistic','fndn90_realistic']
title_exp = ['CTRL','ANTH','50% N Red.','50% N Red.\n50% Recy.','50% N Red.\n90% Recy.','85% N Red.','85% N Red.\n50% Recy.','85% N Red.\n90% Recy.']

# 2 year average
# average y: 450-460, x all points
avg2year_path = '/data/project6/minnaho/opc_scenarios/avg2years/'
avg2year_files = glob.glob(avg2year_path+'*')

axfont=16

# anth vs cntrl
fig,ax = plt.subplots(1,2,figsize=[10,7])
cntrl_nc = Dataset(avg2year_path+'avg2years_'+exp[0]+'.nc','r')
fulll_nc = Dataset(avg2year_path+'avg2years_'+exp[1]+'.nc','r')

z_cntrl = z_r = rd.get_zr_zw_tind(cntrl_nc,grid_nc,0,[450,460,0,cntrl_nc['temp'].shape[3]])[0]
z_fulll = z_r = rd.get_zr_zw_tind(fulll_nc,grid_nc,0,[450,460,0,fulll_nc['temp'].shape[3]])[0]

z_cntrl[z_cntrl>1E10] = np.nan
z_cntrl_plt = np.nanmean(z_cntrl,axis=1)

z_fulll[z_fulll>1E10] = np.nan
z_fulll_plt = np.nanmean(z_fulll,axis=1)

if var_nc != 'biomass':
    data_cntrl = np.squeeze(cntrl_nc[var_nc])
    data_cntrl[data_cntrl>1E10] = np.nan
    slice_cntrl = np.nanmean(data_cntrl[:,450:460,:],axis=1)
    ax.flat[0].pcolor((np.arange(data_cntrl1.shape[2])-575)*-1*xplt,z_cntrl_plt,slice_cntrl)
    ax.flat[0].set_ylim(bottom=-100,top=0)
    ax.flat[0].set_xlim([0,95])
    ax.flat[0].invert_xaxis()

    data_fulll = np.squeeze(fulll_nc[var_nc])
    data_fulll[data_fulll>1E10] = np.nan
    slice_fulll = np.nanmean(data_fulll[:,450:460,:],axis=1)
    ax.flat[1].pcolor((np.arange(data_fulll1.shape[2])-575)*-1*xplt,z_fulll_plt,slice_fulll)
    ax.flat[1].set_ylim(bottom=-100,top=0)
    ax.flat[1].set_xlim([0,95])
    ax.flat[1].invert_xaxis()

if var_nc == 'biomass':
    data_cntrl1 = np.squeeze(cntrl_nc['DIATC'])
    data_cntrl2 = np.squeeze(cntrl_nc['DIAZC'])
    data_cntrl3 = np.squeeze(cntrl_nc['SPC'])
    data_cntrl1[data_cntrl1>1E10] = np.nan
    data_cntrl2[data_cntrl2>1E10] = np.nan
    data_cntrl3[data_cntrl3>1E10] = np.nan
    slice_cntrl = np.nanmean(data_cntrl1[:,450:460,:],axis=1)+np.nanmean(data_cntrl2[:,450:460,:],axis=1)+np.nanmean(data_cntrl3[:,450:460,:],axis=1)
    ax.flat[0].pcolor((np.arange(data_cntrl1.shape[2])-575)*-1*xplt,z_cntrl_plt,slice_cntrl)
    ax.flat[0].set_ylim(bottom=-100,top=0)
    ax.flat[0].set_xlim([0,95])
    ax.flat[0].invert_xaxis()

    data_fulll1 = np.squeeze(fulll_nc['DIATC'])
    data_fulll2 = np.squeeze(fulll_nc['DIAZC'])
    data_fulll3 = np.squeeze(fulll_nc['SPC'])
    data_fulll1[data_fulll1>1E10] = np.nan
    data_fulll2[data_fulll2>1E10] = np.nan
    data_fulll3[data_fulll3>1E10] = np.nan
    slice_fulll = np.nanmean(data_fulll1[:,450:460,:],axis=1)+np.nanmean(data_fulll2[:,450:460,:],axis=1)+np.nanmean(data_fulll3[:,450:460,:],axis=1)
    pplot1 = ax.flat[1].pcolor((np.arange(data_fulll1.shape[2])-575)*-1*xplt,z_fulll_plt,slice_fulll-slice_cntrl,cmap=cmap1)
    ax.flat[1].set_ylim(bottom=-100,top=0)
    ax.flat[1].set_xlim([0,95])
    ax.flat[1].invert_xaxis()
    p0 = ax.flat[1].get_position().get_points().flatten()
    p1 = ax.flat[1].get_position().get_points().flatten()
    cb_ax = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])
    cb1 = fig.colorbar(pplot1,cax=cb_ax,orientation='vertical')
    cb1.set_label('$\Delta$ biomass mmol m$^{-3}$',fontsize=axfont)



'''
fig,ax = plt.subplots(2,3,figsize=[18,10])
for e_i in range(2,len(exp)):
    avg2year_nc = Dataset(avg2year_path+'avg2years_'+exp[e_i]+'.nc','r')
    # calculate depths
    z_r = rd.get_zr_zw_tind(avg2year_nc,grid_nc,0,[450,460,0,avg2year_nc['temp'].shape[3]])[0]
    z_r[z_r>1E10] = np.nan
    z_plt = np.nanmean(z_r,axis=1)
    # calculate mean over y=450-460
    if var_nc != 'biomass':
        datanc = np.squeeze(avg2year_nc[var_nc])
        datanc[datanc>1E10] = np.nan
        ncslice = np.nanmean(datanc[:,450:460,:],axis=1)
        ax.flat[e_i-2].pcolor(np.arange(datanc.shape[2]),z_plt,ncslice)
    if var_nc == 'biomass':
        datanc1 = np.squeeze(avg2year_nc['DIATC'])
        datanc2 = np.squeeze(avg2year_nc['DIAZC'])
        datanc3 = np.squeeze(avg2year_nc['SPC'])
        datanc1[datanc1>1E10] = np.nan
        datanc2[datanc2>1E10] = np.nan
        datanc3[datanc3>1E10] = np.nan
        ncslice = np.nanmean(datanc1[:,450:460,:],axis=1)+np.nanmean(datanc2[:,450:460,:],axis=1)+np.nanmean(datanc3[:,450:460,:],axis=1)
        ax.flat[e_i-2].pcolor(np.arange(datanc1.shape[2]),z_plt,ncslice)
        ax.flat[e_i-2].set_ylim(bottom=-100,top=0)
'''    
