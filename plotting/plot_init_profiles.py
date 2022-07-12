import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import subprocess
import l2grid 
import ROMS_depths as rd
import pyroms

# check/convert to actual depths
# not sure how to handle initial 

outpath = '/data/project6/minnaho/opc_scenarios/plume_shape/'

varstr = 'NH4'
dist = 1
dst = 1
den = 90
stp = 1

# avg and std of monthly files
'''
var = 'DIATC,DIAZC,SPC,zeta'
exp = 'pndn90'
romspath = '/data/project6/ROMS/L2SCB_OPC/'+exp+'/monthly/'

print('concat '+exp)
subprocess.call('ncrcat -v '+var+' '+romspath+'l2_scb_avg.* '+outpath+'concat_'+var+'_'+exp+'.nc',shell=True)
print('avg '+exp)
subprocess.call('ncwa -a time '+outpath+'concat_'+var+'_'+exp+'.nc'+' '+outpath+'avg_'+var+'_'+exp+'.nc',shell=True)
print('anom '+exp)
subprocess.call('ncbo -v '+var+' '+outpath+'concat_'+var+'_'+exp+'.nc'+' '+outpath+'avg_'+var+'_'+exp+'.nc '+outpath+'anom_'+var+'_'+exp+'.nc',shell=True)
print('std '+exp)
subprocess.call('ncra -y rmssdn '+outpath+'anom_'+var+'_'+exp+'.nc'+' '+outpath+'std_'+var+'_'+exp+'.nc',shell=True)
'''

# get grid
grd = pyroms.grid.get_ROMS_grid('L2')

gridnc = l2grid.grid_nc


# read in major POTW radii
outf = ['hyp','jwp','ocs','plw']

# range of depths
dps = np.arange(dst,den+stp,stp)


# plot profiles of plume
# PNDN 
if varstr == 'NH4':
    pndnon_avg_read = np.squeeze(Dataset(outpath+'avg_NH4_PNDN_only.nc','r').variables['NH4'])
    pndnon_std_read = np.squeeze(Dataset(outpath+'std_NH4_PNDN_only.nc','r').variables['NH4'])
    # pndn50 
    pndn50_avg_read = np.squeeze(Dataset(outpath+'avg_NH4_pndn50.nc','r').variables['NH4'])
    pndn50_std_read = np.squeeze(Dataset(outpath+'std_NH4_pndn50.nc','r').variables['NH4'])
    # pndn90
    pndn90_avg_read = np.squeeze(Dataset(outpath+'avg_NH4_pndn90.nc','r').variables['NH4'])
    pndn90_std_read = np.squeeze(Dataset(outpath+'std_NH4_pndn90.nc','r').variables['NH4'])
    hline = 3


if varstr == 'biomass':
    pndnon_avg_read = np.squeeze(Dataset(outpath+'avg_DIATC_DIAZC_SPC_zeta_PNDN_only.nc','r').variables['DIATC'])+np.squeeze(Dataset(outpath+'avg_DIATC_DIAZC_SPC_zeta_PNDN_only.nc','r').variables['DIAZC'])+np.squeeze(Dataset(outpath+'avg_DIATC_DIAZC_SPC_zeta_PNDN_only.nc','r').variables['SPC'])
    pndnon_std_read = np.squeeze(Dataset(outpath+'avg_DIATC_DIAZC_SPC_zeta_PNDN_only.nc','r').variables['DIATC'])+np.squeeze(Dataset(outpath+'avg_DIATC_DIAZC_SPC_zeta_PNDN_only.nc','r').variables['DIAZC'])+np.squeeze(Dataset(outpath+'avg_DIATC_DIAZC_SPC_zeta_PNDN_only.nc','r').variables['SPC'])
    pndn50_avg_read = np.squeeze(Dataset(outpath+'avg_DIATC_DIAZC_SPC_zeta_pndn50.nc','r').variables['DIATC'])+np.squeeze(Dataset(outpath+'avg_DIATC_DIAZC_SPC_zeta_pndn50.nc','r').variables['DIAZC'])+np.squeeze(Dataset(outpath+'avg_DIATC_DIAZC_SPC_zeta_pndn50.nc','r').variables['SPC'])
    pndn50_std_read = np.squeeze(Dataset(outpath+'avg_DIATC_DIAZC_SPC_zeta_pndn50.nc','r').variables['DIATC'])+np.squeeze(Dataset(outpath+'avg_DIATC_DIAZC_SPC_zeta_pndn50.nc','r').variables['DIAZC'])+np.squeeze(Dataset(outpath+'avg_DIATC_DIAZC_SPC_zeta_pndn50.nc','r').variables['SPC'])
    pndn90_avg_read = np.squeeze(Dataset(outpath+'avg_DIATC_DIAZC_SPC_zeta_pndn90.nc','r').variables['DIATC'])+np.squeeze(Dataset(outpath+'avg_DIATC_DIAZC_SPC_zeta_pndn90.nc','r').variables['DIAZC'])+np.squeeze(Dataset(outpath+'avg_DIATC_DIAZC_SPC_zeta_pndn90.nc','r').variables['SPC'])
    pndn90_std_read = np.squeeze(Dataset(outpath+'avg_DIATC_DIAZC_SPC_zeta_pndn90.nc','r').variables['DIATC'])+np.squeeze(Dataset(outpath+'avg_DIATC_DIAZC_SPC_zeta_pndn90.nc','r').variables['DIAZC'])+np.squeeze(Dataset(outpath+'avg_DIATC_DIAZC_SPC_zeta_pndn90.nc','r').variables['SPC'])
    hline = 10

pndnon_avg_read[pndnon_avg_read>1E10] = np.nan
pndnon_std_read[pndnon_std_read>1E10] = np.nan

pndn50_avg_read[pndn50_avg_read>1E10] = np.nan
pndn50_std_read[pndn50_std_read>1E10] = np.nan

pndn90_avg_read[pndn90_avg_read>1E10] = np.nan
pndn90_std_read[pndn90_std_read>1E10] = np.nan

# concatenate into one file to loop
rfiles_avg = np.ones((3,pndnon_avg_read.shape[0],pndnon_avg_read.shape[1],pndnon_avg_read.shape[2]))*np.nan

rfiles_std = np.ones((3,pndnon_avg_read.shape[0],pndnon_avg_read.shape[1],pndnon_avg_read.shape[2]))*np.nan

rfiles_avg[0,:,:,:] = pndnon_avg_read[:,:,:]
rfiles_avg[1,:,:,:] = pndn50_avg_read[:,:,:]
rfiles_avg[2,:,:,:] = pndn90_avg_read[:,:,:]

rfiles_std[0,:,:,:] = pndnon_std_read[:,:,:]
rfiles_std[1,:,:,:] = pndn50_std_read[:,:,:]
rfiles_std[2,:,:,:] = pndn90_std_read[:,:,:]

prof_avg = np.ones((rfiles_avg.shape[0],4,dps.shape[0]))*np.nan
prof_std = np.ones((rfiles_avg.shape[0],4,dps.shape[0]))*np.nan

# initial profile
init = '/data/project1/minnaho/psource/wastewater_scenarios/roms_psource_PNDN_only_realistic.nc'
htp_st = 0
jwp_st = 28
ocs_st = 56
plw_st = 70
plw_en = 96

qshp = np.squeeze(Dataset(init,'r').variables['Qshape'])
isrc = np.squeeze(Dataset(init,'r').variables['Isrc'])
jsrc = np.squeeze(Dataset(init,'r').variables['Jsrc'])
nval = np.squeeze(Dataset(init,'r').variables['NH4'])

htpq = qshp[:,htp_st:jwp_st]
jwpq = qshp[:,jwp_st:ocs_st]
ocsq = qshp[:,ocs_st:plw_st]
plwq = qshp[:,plw_st:plw_en]

htpi = isrc[htp_st:jwp_st].astype(int)
jwpi = isrc[jwp_st:ocs_st].astype(int)
ocsi = isrc[ocs_st:plw_st].astype(int)
plwi = isrc[plw_st:plw_en].astype(int)
                          
htpj = jsrc[htp_st:jwp_st].astype(int)
jwpj = jsrc[jwp_st:ocs_st].astype(int)
ocsj = jsrc[ocs_st:plw_st].astype(int)
plwj = jsrc[plw_st:plw_en].astype(int)

htpn = nval[htp_st]
jwpn = nval[jwp_st]
ocsn = nval[ocs_st]
plwn = nval[plw_st]

pndnon_avg_zeta = np.squeeze(Dataset(outpath+'avg_NH4_zeta_PNDN_only.nc','r').variables['zeta'])


expname = 'PNDN_only'
fnc = Dataset(outpath+'avg_NH4_zeta_PNDN_only.nc','r')
zdepths = rd.get_zs3d(fnc,gridnc,dim_bounds=[0,pndnon_avg_zeta.shape[0],0,pndnon_avg_zeta.shape[1]],varname='NH4')
zdepths[zdepths>1E10] = np.nan

htpz = zdepths[:,htpj,htpi]
jwpz = zdepths[:,jwpj,jwpi]
ocsz = zdepths[:,ocsj,ocsi]
plwz = zdepths[:,plwj,plwi]

plt.ion()

figw = 8
figh = 10

hline = 50

fig,ax = plt.subplots(1,1,figsize=[figw,figh])
for p_i in range(htpz.shape[1]):
    ax.plot(htpq[:,p_i]*htpn[0],htpz.T[p_i]*-1)
    ax.plot(np.arange(hline),np.ones((len(np.arange(hline))))*htpz.T[p_i][np.where(htpq[:,p_i]*htpn[0]==np.nanmax(htpq[:,p_i]*htpn[0]))[0][0]]*-1,color='blue',linestyle=':')
    ax.set_title(outf[0]+' initial')
    ax.set_ylabel('Depth (m)')
    ax.set_xlabel(varstr+' (mmol/m3)')
    ax.set_ylim([-3,den+3])
    ax.invert_yaxis()
    ax.legend(loc='best')
    fig.savefig('./figs/profiles/initial_'+outf[0]+'_'+varstr+'_'+expname+'.png',bbox_inches='tight')

fig,ax = plt.subplots(1,1,figsize=[figw,figh])
for p_i in range(jwpz.shape[1]):
    ax.plot(jwpq[:,p_i]*jwpn[0],jwpz.T[p_i]*-1)
    ax.plot(np.arange(hline),np.ones((len(np.arange(hline))))*jwpz.T[p_i][np.where(jwpq[:,p_i]*jwpn[0]==np.nanmax(jwpq[:,p_i]*jwpn[0]))[0][0]]*-1,color='blue',linestyle=':')
    ax.set_title(outf[1]+' initial')
    ax.set_ylabel('Depth (m)')
    ax.set_xlabel(varstr+' (mmol/m3)')
    ax.set_ylim([-3,den+3])
    ax.invert_yaxis()
    ax.legend(loc='best')
    fig.savefig('./figs/profiles/initial_'+outf[1]+'_'+varstr+'_'+expname+'.png',bbox_inches='tight')


fig,ax = plt.subplots(1,1,figsize=[figw,figh])
for p_i in range(ocsz.shape[1]):
    ax.plot(ocsq[:,p_i]*ocsn[0],ocsz.T[p_i]*-1)
    ax.plot(np.arange(hline),np.ones((len(np.arange(hline))))*ocsz.T[p_i][np.where(ocsq[:,p_i]*ocsn[0]==np.nanmax(ocsq[:,p_i]*ocsn[0]))[0][0]]*-1,color='blue',linestyle=':')
    ax.set_title(outf[2]+' initial')
    ax.set_ylabel('Depth (m)')
    ax.set_xlabel(varstr+' (mmol/m3)')
    ax.set_ylim([-3,den+3])
    ax.invert_yaxis()
    ax.legend(loc='best')
    fig.savefig('./figs/profiles/initial_'+outf[2]+'_'+varstr+'_'+expname+'.png',bbox_inches='tight')

fig,ax = plt.subplots(1,1,figsize=[figw,figh])
for p_i in range(plwz.shape[1]):
    ax.plot(plwq[:,p_i]*plwn[0],plwz.T[p_i]*-1)
    ax.plot(np.arange(hline),np.ones((len(np.arange(hline))))*plwz.T[p_i][np.where(plwq[:,p_i]*plwn[0]==np.nanmax(plwq[:,p_i]*plwn[0]))[0][0]]*-1,color='blue',linestyle=':')
    ax.set_title(outf[3]+' initial')
    ax.set_ylabel('Depth (m)')
    ax.set_xlabel(varstr+' (mmol/m3)')
    ax.set_ylim([-3,den+3])
    ax.invert_yaxis()
    ax.legend(loc='best')
    fig.savefig('./figs/profiles/initial_'+outf[3]+'_'+varstr+'_'+expname+'.png',bbox_inches='tight')

