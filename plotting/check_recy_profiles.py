import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import subprocess

# check/convert to actual depths
# not sure how to handle initial 

outpath = '/data/project6/minnaho/opc_scenarios/plume_shape/'

# avg and std of monthly files
'''
var = 'NH4'
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

# read in major POTW radii
hypr = np.squeeze(Dataset(outpath+'hyp5km.nc','r').variables['mask'])
jwpr = np.squeeze(Dataset(outpath+'jwp5km.nc','r').variables['mask'])
ocsr = np.squeeze(Dataset(outpath+'ocs5km.nc','r').variables['mask'])
plwr = np.squeeze(Dataset(outpath+'plw5km.nc','r').variables['mask'])

hypr[hypr==0] = np.nan
jwpr[jwpr==0] = np.nan
ocsr[ocsr==0] = np.nan
plwr[plwr==0] = np.nan

# combine into one
#hypr[np.where(jwpr==1)[0],np.where(jwpr==1)[1]] = 1
#hypr[np.where(ocsr==1)[0],np.where(ocsr==1)[1]] = 1
#hypr[np.where(plwr==1)[0],np.where(plwr==1)[1]] = 1

# choose mask
maskch = hypr

# plot profiles of plume
# PNDN 
pndnon_avg_read = np.squeeze(Dataset(outpath+'avg_NH4_PNDN_only.nc','r').variables['NH4'])*maskch
pndnon_std_read = np.squeeze(Dataset(outpath+'std_NH4_PNDN_only.nc','r').variables['NH4'])*maskch
# pndn50 
pndn50_avg_read = np.squeeze(Dataset(outpath+'avg_NH4_pndn50.nc','r').variables['NH4'])*maskch
pndn50_std_read = np.squeeze(Dataset(outpath+'std_NH4_pndn50.nc','r').variables['NH4'])*maskch
# pndn90
pndn90_avg_read = np.squeeze(Dataset(outpath+'avg_NH4_pndn90.nc','r').variables['NH4'])*maskch
pndn90_std_read = np.squeeze(Dataset(outpath+'std_NH4_pndn90.nc','r').variables['NH4'])*maskch

pndnon_avg_read[pndnon_avg_read>1E10] = np.nan
pndnon_std_read[pndnon_std_read>1E10] = np.nan

pndn50_avg_read[pndn50_avg_read>1E10] = np.nan
pndn50_std_read[pndn50_std_read>1E10] = np.nan

pndn90_avg_read[pndn90_avg_read>1E10] = np.nan
pndn90_std_read[pndn90_std_read>1E10] = np.nan

# average over mask

pndnon_avg = np.squeeze(np.apply_over_axes(np.nanmean,pndnon_avg_read,(1,2)))
pndnon_std = np.squeeze(np.apply_over_axes(np.nanmean,pndnon_std_read,(1,2)))
# pndn50                                               
pndn50_avg = np.squeeze(np.apply_over_axes(np.nanmean,pndn50_avg_read,(1,2)))
pndn50_std = np.squeeze(np.apply_over_axes(np.nanmean,pndn50_std_read,(1,2)))
# pndn90                                               
pndn90_avg = np.squeeze(np.apply_over_axes(np.nanmean,pndn90_avg_read,(1,2)))
pndn90_std = np.squeeze(np.apply_over_axes(np.nanmean,pndn90_std_read,(1,2)))

# initial profile
init = '/data/project1/minnaho/psource/wastewater_scenarios/roms_psource_PNDN_only_realistic.nc'
htp_st = 0
jwp_st = 28
ocs_st = 56
plw_st = 70

qshp = np.squeeze(Dataset(init,'r').variables['Qshape'])
htpq = qshp[:,htp_st]
jwpq = qshp[:,jwp_st]
ocsq = qshp[:,ocs_st]
plwq = qshp[:,plw_st]

srho = np.arange(qshp.shape[0])

figw = 5
figh = 10

plt.ion()

fig,ax = plt.subplots(1,1,figsize=[figw,figh])
ax.plot(htpq,srho,label='initial',color='black')
ax.plot(pndnon_avg,srho,label='PNDN',color='blue')
ax.plot(pndn50_avg,srho,label='PNDN50',color='blue',linestyle='--')
ax.plot(pndn90_avg,srho,label='PNDN90',color='lightblue',linestyle=':')

