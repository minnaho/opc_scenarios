################################################
# plot cross section of freshwater/nutrient/control/full 
# at major pipes
################################################
import sys
import os
sys.path.append('/data/project3/minnaho/global/')
import l2grid
import numpy as np
from netCDF4 import Dataset,num2date
import glob as glob
import ROMS_depths as depths
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import cmocean
import datetime as datetime
import calendar
import seawater as sw

# plot fresh vs nutrients vs control vs full
#plt.ion()

savepath = './figs/cs/'
loc = 'OCSan'

# roms var
var_nc = 'N2' # buoyancy frequency
#var_nc = 'biomass'
#cblabel = 'Density (kg m$^{-3}$)'
#cblabel = 'Temperature C'
#cblabel = 'NO3 (mmol m$^{-3}$)'
#cblabel = 'NH4 (mmol m$^{-3}$)'
#cblabel = 'mmol C m$^{-3}$'
cblabel = 'N$^2$ s$^{-1}$'
#cblabel = 'm s$^{-1}$'
#cblabel = 'salt (PSU)'

# path of outputs
cntrlpath = '/data/project6/ROMS/L2SCB_1997_2000/monthly/'
fulllpath = '/data/project6/ROMS/L2SCB_AP/monthly/'
lo1617path = '/data/project6/ROMS/L2SCB_OPC/loads1617/monthly/'
pndnonpath = '/data/project6/ROMS/L2SCB_OPC/PNDN_only/monthly/'
fndnonpath = '/data/project6/ROMS/L2SCB_OPC/FNDN_only/monthly/'
pndn50path = '/data/project6/ROMS/L2SCB_OPC/pndn50/monthly/'
pndn90path = '/data/project6/ROMS/L2SCB_OPC/pndn90/monthly/'
fndn50path = '/data/project6/ROMS/L2SCB_OPC/fndn50/monthly/'
fndn90path = '/data/project6/ROMS/L2SCB_OPC/fndn90/monthly/'

start_year = 1999
end_year = 2000

# between 1 and 12
start_month = 9
end_month = 9


if var_nc == 'rho':
    rho0 = 1027.4
    clines = [1023.5,1024,1024.5,1025,1025.5,1026,1026.5]
if var_nc == 'NH4' or var_nc == 'NO3':
    clines = [0.5,1,5,10,13]
if var_nc == 'salt':
    clines = [32.75,33,33.25,33.5,33.75,34,34.25]
if var_nc == 'temp':
    clines = [9,10,11,12,13,14,15,16,17,18,19,20]
if var_nc == 'biomass':
    clines = [0.5,1,3,5,10,15,20,25]

#savename = 'cs_'+var_nc+'_Y'+str(year)+'M'+'%02d'%month+'.png'

# roms file name
#ncfile = 'l2_scb_avg.Y'+str(year)+'M'+'%02d'%month+'.nc'

ncfile = []
savename = []
for y in range(start_year,end_year+1):
    # if we are on the first year, starts at s_m
    if y == start_year:
        s_m = start_month
    else:
        s_m = 1
    # if we are on the last year, end at e_m
    if y == end_year:
        e_m = end_month+1
    else:
        e_m = 13
    for m in range(s_m,e_m):
        year_month = 'Y'+str(y)+'M'+'%02d'%m
        ncfile.append('l2_scb_avg.'+year_month+'.nc')
        savename.append('cs_'+loc+'_'+var_nc+'_Y'+str(y)+'M'+'%02d'%m+'.png')

# outputs
cntrlnc = []
fulllnc = []
lo1617nc = []
pndnonnc = []
fndnonnc = []
pndn50nc = []
pndn90nc = []
fndn50nc = []
fndn90nc = []

for n_i in range(len(ncfile)):
    cntrlnc.append(cntrlpath+ncfile[n_i])
    fulllnc.append(fulllpath+ncfile[n_i])
    lo1617nc.append(lo1617path+ncfile[n_i])
    pndnonnc.append(pndnonpath+ncfile[n_i])
    fndnonnc.append(fndnonpath+ncfile[n_i])
    pndn50nc.append(pndn50path+ncfile[n_i])
    pndn90nc.append(pndn90path+ncfile[n_i])
    fndn50nc.append(fndn50path+ncfile[n_i])
    fndn90nc.append(fndn90path+ncfile[n_i])



# grid path
grid_nc = l2grid.grid_nc
lat_nc = l2grid.lat_nc
lon_nc = l2grid.lon_nc
h_nc = l2grid.h_nc

if loc == 'HTP':
    # north pipe
    lat_site = 33.920625
    lon_site = -118.529712
    ind_st = -118.6
    ind_en = -118.46
    
    dp_st = -80
    dp_en = 0

    dpipe = -60
    
    #p_loc = 551 # location to mark pipe

if loc == 'JWPCP':
    # 65% diffuser
    #lat_site = 33.6892
    #lon_site = -118.3167

    # 35% diffuser
    lat_site = 33.700737
    lon_site = -118.341962
    ind_st = -118.43
    ind_en = -118.31
    
    dp_st = -80
    dp_en = 0

    dpipe = -50

# ocsd pipe - take cross section here
if loc == 'OCSan':
    lat_site = 33.57667
    lon_site = -118.01
    #ind_st_p = 500 # remove all offshore
    #ind_en_p = 571 # remove land
    ind_st = -118.1
    ind_en = -117.963
    
    dp_st = -70
    dp_en = 0
    
    p_loc = 551 # location to mark pipe

if loc == 'PLWTP':
    # north pipe
    lat_site = 32.671671
    lon_site = -117.325556 
    ind_st = -117.4
    ind_en = -117.25
    
    dp_st = -110
    dp_en = 0

    dpipe = -95


# calculate i and j
min_1D = np.abs( (lat_nc - lat_site)**2 + (lon_nc - lon_site)**2)
y_site, x_site = np.unravel_index(min_1D.argmin(), min_1D.shape)

# get longitudinal slice
lon_slice = lon_nc[y_site,:]
h_slice = h_nc[y_site,:]

# reshape to match z_r
# 60 by 602 in L2
lon_slice_l = list(lon_slice)*Dataset(freshnc[0],'r').variables['temp'].shape[1]
lon_reshape = np.array(lon_slice_l).reshape(Dataset(freshnc[0],'r').variables['temp'].shape[1],lon_nc.shape[1]) 

# bounds to draw contours
ind_st_p = np.nanmin(np.unique(np.where((lon_reshape[:,:]>ind_st)&(lon_reshape<ind_en))[1]))
ind_en_p = np.nanmax(np.unique(np.where((lon_reshape[:,:]>ind_st)&(lon_reshape<ind_en))[1]))

figw = 14
figh = 7.5

if var_nc == 'w':
    c_map = cmocean.cm.balance
else:
    c_map = cmocean.cm.dense

# size of marker for pipe location
psize = 300

axis_tick_size = 14

for n_i in range(len(ncfile)):
    fig,ax = plt.subplots(2,4,figsize=[figw,figh])

    # get depths at this y_site slice for each scenario
    z_r_cntrl = depths.get_zr_zw_tind(Dataset(cntrlnc[n_i],'r'),grid_nc,0,[y_site-1,y_site+1,0,Dataset(cntrlnc[n_i],'r').variables['temp'].shape[3]])[0][:,1,:]
    #z_r_fulll = depths.get_zr_zw_tind(Dataset(fulllnc[n_i],'r'),grid_nc,0,[y_site-1,y_site+1,0,Dataset(fulllnc[n_i],'r').variables['temp'].shape[3]])[0][:,1,:]
    z_r_lo1617 = depths.get_zr_zw_tind(Dataset(lo1617nc[n_i],'r'),grid_nc,0,[y_site-1,y_site+1,0,Dataset(lo1617nc[n_i],'r').variables['temp'].shape[3]])[0][:,1,:]
    z_r_pndnon = depths.get_zr_zw_tind(Dataset(pndnonnc[n_i],'r'),grid_nc,0,[y_site-1,y_site+1,0,Dataset(pndnonnc[n_i],'r').variables['temp'].shape[3]])[0][:,1,:]
    z_r_fndnon = depths.get_zr_zw_tind(Dataset(fndnonnc[n_i],'r'),grid_nc,0,[y_site-1,y_site+1,0,Dataset(fndnonnc[n_i],'r').variables['temp'].shape[3]])[0][:,1,:]
    z_r_pndn50 = depths.get_zr_zw_tind(Dataset(pndn50nc[n_i],'r'),grid_nc,0,[y_site-1,y_site+1,0,Dataset(pndn50nc[n_i],'r').variables['temp'].shape[3]])[0][:,1,:]
    z_r_pndn90 = depths.get_zr_zw_tind(Dataset(pndn90nc[n_i],'r'),grid_nc,0,[y_site-1,y_site+1,0,Dataset(pndn90nc[n_i],'r').variables['temp'].shape[3]])[0][:,1,:]
    z_r_fndn50 = depths.get_zr_zw_tind(Dataset(fndn50nc[n_i],'r'),grid_nc,0,[y_site-1,y_site+1,0,Dataset(fndn50nc[n_i],'r').variables['temp'].shape[3]])[0][:,1,:]
    z_r_fndn90 = depths.get_zr_zw_tind(Dataset(fndn90nc[n_i],'r'),grid_nc,0,[y_site-1,y_site+1,0,Dataset(fndn90nc[n_i],'r').variables['temp'].shape[3]])[0][:,1,:]

    # roms field from output at y_site slice
    if var_nc == 'rho':
        roms_var_cntrl = np.array(Dataset(cntrlnc[n_i],'r').variables[var_nc][0,:,y_site,:])+rho0
        #roms_var_fulll = np.array(Dataset(fulllnc[n_i],'r').variables[var_nc][0,:,y_site,:])+rho0
        roms_var_lo1617 = np.array(Dataset(lo1617nc[n_i],'r').variables[var_nc][0,:,y_site,:])+rho0
        roms_var_pndnon = np.array(Dataset(pndnonnc[n_i],'r').variables[var_nc][0,:,y_site,:])+rho0
        roms_var_fndnon = np.array(Dataset(fndnonnc[n_i],'r').variables[var_nc][0,:,y_site,:])+rho0
        roms_var_pndn50 = np.array(Dataset(pndn50nc[n_i],'r').variables[var_nc][0,:,y_site,:])+rho0
        roms_var_pndn90 = np.array(Dataset(pndn90nc[n_i],'r').variables[var_nc][0,:,y_site,:])+rho0
        roms_var_fndn50 = np.array(Dataset(fndn50nc[n_i],'r').variables[var_nc][0,:,y_site,:])+rho0
        roms_var_fndn90 = np.array(Dataset(fndn90nc[n_i],'r').variables[var_nc][0,:,y_site,:])+rho0
    if var_nc == 'N2':
        # calculate N2 brunt vaisalla buoyancy frequency 
        # salt
        cntrlsal = np.array(Dataset(cntrlnc[n_i],'r').variables['salt'][0,:,y_site,:])
        #fulllsal = np.array(Dataset(fulllnc[n_i],'r').variables['salt'][0,:,y_site,:])
        lo1617sal = np.array(Dataset(lo1617nc[n_i],'r').variables['salt'][0,:,y_site,:])
        pndnonsal = np.array(Dataset(pndnonnc[n_i],'r').variables['salt'][0,:,y_site,:])
        fndnonsal = np.array(Dataset(fndnonnc[n_i],'r').variables['salt'][0,:,y_site,:])
        pndn50sal = np.array(Dataset(pndn50nc[n_i],'r').variables['salt'][0,:,y_site,:])
        pndn90sal = np.array(Dataset(pndn90nc[n_i],'r').variables['salt'][0,:,y_site,:])
        fndn50sal = np.array(Dataset(fndn50nc[n_i],'r').variables['salt'][0,:,y_site,:])
        fndn90sal = np.array(Dataset(fndn90nc[n_i],'r').variables['salt'][0,:,y_site,:])

        # temp
        cntrltem = np.array(Dataset(cntrlnc[n_i],'r').variables['temp'][0,:,y_site,:])
        #fullltem = np.array(Dataset(fulllnc[n_i],'r').variables['temp'][0,:,y_site,:])
        lo1617tem = np.array(Dataset(lo1617nc[n_i],'r').variables['temp'][0,:,y_site,:])
        pndnontem = np.array(Dataset(pndnonnc[n_i],'r').variables['temp'][0,:,y_site,:])
        fndnontem = np.array(Dataset(fndnonnc[n_i],'r').variables['temp'][0,:,y_site,:])
        pndn50tem = np.array(Dataset(pndn50nc[n_i],'r').variables['temp'][0,:,y_site,:])
        pndn90tem = np.array(Dataset(pndn90nc[n_i],'r').variables['temp'][0,:,y_site,:])
        fndn50tem = np.array(Dataset(fndn50nc[n_i],'r').variables['temp'][0,:,y_site,:])
        fndn90tem = np.array(Dataset(fndn90nc[n_i],'r').variables['temp'][0,:,y_site,:])

        # pressure
        cntrlpre = sw.pres(z_r_cntrl*-1,lat_site)
        #fulllpre = sw.pres(z_r_fulll*-1,lat_site)
        lo1617pre = sw.pres(z_r_lo1617*-1,lat_site)
        pndnonpre = sw.pres(z_r_pndnon*-1,lat_site)
        fndnonpre = sw.pres(z_r_fndnon*-1,lat_site)
        pndn50pre = sw.pres(z_r_pndn50*-1,lat_site)
        pndn90pre = sw.pres(z_r_pndn90*-1,lat_site)
        fndn50pre = sw.pres(z_r_fndn50*-1,lat_site)
        fndn90pre = sw.pres(z_r_fndn90*-1,lat_site)

        roms_var_cntrl = sw.bfrq(cntrlsal,cntrltem,cntrlpre,lat_site)[0]
        #roms_var_fulll = sw.bfrq(fulllsal,fullltem,fulllpre,lat_site)[0]
        roms_var_lo1617 = sw.bfrq(lo1617sal,lo1617tem,lo1617pre,lat_site)[0]
        roms_var_pndnon = sw.bfrq(pndnonsal,pndnontem,pndnonpre,lat_site)[0]
        roms_var_fndnon = sw.bfrq(fndnonsal,fndnontem,fndnonpre,lat_site)[0]
        roms_var_pndn50 = sw.bfrq(pndn50sal,pndn50tem,pndn50pre,lat_site)[0]
        roms_var_pndn90 = sw.bfrq(pndn90sal,pndn90tem,pndn90pre,lat_site)[0]
        roms_var_fndn50 = sw.bfrq(fndn50sal,fndn50tem,fndn50pre,lat_site)[0]
        roms_var_fndn90 = sw.bfrq(fndn90sal,fndn90tem,fndn90pre,lat_site)[0]

    if var_nc == 'DIN':
        roms_var_cntrl_nh4 = np.array(Dataset(cntrlnc[n_i],'r').variables['NH4'][0,:,y_site,:])
        roms_var_cntrl_nh4[roms_var_cntrl_nh4>1E10] = 0
        roms_var_cntrl_no3 = np.array(Dataset(cntrlnc[n_i],'r').variables['NO3'][0,:,y_site,:])
        roms_var_cntrl_no3[roms_var_cntrl_no3>1E10] = 0

        #roms_var_fulll_nh4 = np.array(Dataset(fulllnc[n_i],'r').variables['NH4'][0,:,y_site,:])
        #roms_var_fulll_nh4[roms_var_fulll_nh4>1E10] = 0
        #roms_var_fulll_no3 = np.array(Dataset(fulllnc[n_i],'r').variables['NO3'][0,:,y_site,:])
        #roms_var_fulll_no3[roms_var_fulll_no3>1E10] = 0

        roms_var_lo1617_nh4 = np.array(Dataset(lo1617nc[n_i],'r').variables['NH4'][0,:,y_site,:])
        roms_var_lo1617_nh4[roms_var_lo1617_nh4>1E10] = 0
        roms_var_lo1617_no3 = np.array(Dataset(lo1617nc[n_i],'r').variables['NO3'][0,:,y_site,:])
        roms_var_lo1617_no3[roms_var_lo1617_no3>1E10] = 0

        roms_var_pndnon_nh4 = np.array(Dataset(pndnonnc[n_i],'r').variables['NH4'][0,:,y_site,:])
        roms_var_pndnon_nh4[roms_var_pndnon_nh4>1E10] = 0
        roms_var_pndnon_no3 = np.array(Dataset(pndnonnc[n_i],'r').variables['NO3'][0,:,y_site,:])
        roms_var_pndnon_no3[roms_var_pndnon_no3>1E10] = 0

        roms_var_fndnon_nh4 = np.array(Dataset(fndnonnc[n_i],'r').variables['NH4'][0,:,y_site,:])
        roms_var_fndnon_nh4[roms_var_fndnon_nh4>1E10] = 0
        roms_var_fndnon_no3 = np.array(Dataset(fndnonnc[n_i],'r').variables['NO3'][0,:,y_site,:])
        roms_var_fndnon_no3[roms_var_fndnon_no3>1E10] = 0


        roms_var_pndn50_nh4 = np.array(Dataset(pndn50nc[n_i],'r').variables['NH4'][0,:,y_site,:])
        roms_var_pndn50_nh4[roms_var_pndn50_nh4>1E10] = 0
        roms_var_pndn50_no3 = np.array(Dataset(pndn50nc[n_i],'r').variables['NO3'][0,:,y_site,:])
        roms_var_pndn50_no3[roms_var_pndn50_no3>1E10] = 0


        roms_var_pndn90_nh4 = np.array(Dataset(pndn90nc[n_i],'r').variables['NH4'][0,:,y_site,:])
        roms_var_pndn90_nh4[roms_var_pndn90_nh4>1E10] = 0
        roms_var_pndn90_no3 = np.array(Dataset(pndn90nc[n_i],'r').variables['NO3'][0,:,y_site,:])
        roms_var_pndn90_no3[roms_var_pndn90_no3>1E10] = 0


        roms_var_fndn50_nh4 = np.array(Dataset(fndn50nc[n_i],'r').variables['NH4'][0,:,y_site,:])
        roms_var_fndn50_nh4[roms_var_fndn50_nh4>1E10] = 0
        roms_var_fndn50_no3 = np.array(Dataset(fndn50nc[n_i],'r').variables['NO3'][0,:,y_site,:])
        roms_var_fndn50_no3[roms_var_fndn50_no3>1E10] = 0

        roms_var_fndn90_nh4 = np.array(Dataset(fndn90nc[n_i],'r').variables['NH4'][0,:,y_site,:])
        roms_var_fndn90_nh4[roms_var_fndn90_nh4>1E10] = 0
        roms_var_fndn90_no3 = np.array(Dataset(fndn90nc[n_i],'r').variables['NO3'][0,:,y_site,:])
        roms_var_fndn90_no3[roms_var_fndn90_no3>1E10] = 0

        roms_var_cntrl = roms_var_cntrl_dtc+roms_var_cntrl_spc+roms_var_cntrl_dzc
        #roms_var_fulll = roms_var_fulll_dtc+roms_var_fulll_spc+roms_var_fulll_dzc
        roms_var_lo1617 = roms_var_lo1617_dtc+roms_var_lo1617_spc+roms_var_lo1617_dzc
        roms_var_pndnon = roms_var_pndnon_dtc+roms_var_pndnon_spc+roms_var_pndnon_dzc
        roms_var_fndnon = roms_var_fndnon_dtc+roms_var_fndnon_spc+roms_var_fndnon_dzc
        roms_var_pndn50 = roms_var_pndn50_dtc+roms_var_pndn50_spc+roms_var_pndn50_dzc
        roms_var_pndn90 = roms_var_pndn90_dtc+roms_var_pndn90_spc+roms_var_pndn90_dzc
        roms_var_fndn50 = roms_var_fndn50_dtc+roms_var_fndn50_spc+roms_var_fndn50_dzc
        roms_var_fndn90 = roms_var_fndn90_dtc+roms_var_fndn90_spc+roms_var_fndn90_dzc

    if var_nc == 'biomass':
        roms_var_cntrl_dtc = np.array(Dataset(cntrlnc[n_i],'r').variables['DIATC'][0,:,y_site,:])
        roms_var_cntrl_dtc[roms_var_cntrl_dtc>1E10] = 0
        roms_var_cntrl_spc = np.array(Dataset(cntrlnc[n_i],'r').variables['SPC'][0,:,y_site,:])
        roms_var_cntrl_spc[roms_var_cntrl_spc>1E10] = 0
        roms_var_cntrl_dzc = np.array(Dataset(cntrlnc[n_i],'r').variables['DIAZC'][0,:,y_site,:])
        roms_var_cntrl_dzc[roms_var_cntrl_dzc>1E10] = 0

        #roms_var_fulll_dtc = np.array(Dataset(fulllnc[n_i],'r').variables['DIATC'][0,:,y_site,:])
        #roms_var_fulll_dtc[roms_var_fulll_dtc>1E10] = 0
        #roms_var_fulll_spc = np.array(Dataset(fulllnc[n_i],'r').variables['SPC'][0,:,y_site,:])
        #roms_var_fulll_spc[roms_var_fulll_spc>1E10] = 0
        #roms_var_fulll_dzc = np.array(Dataset(fulllnc[n_i],'r').variables['DIAZC'][0,:,y_site,:])
        #roms_var_fulll_dzc[roms_var_fulll_dzc>1E10] = 0

        roms_var_lo1617_dtc = np.array(Dataset(lo1617nc[n_i],'r').variables['DIATC'][0,:,y_site,:])
        roms_var_lo1617_dtc[roms_var_lo1617_dtc>1E10] = 0
        roms_var_lo1617_spc = np.array(Dataset(lo1617nc[n_i],'r').variables['SPC'][0,:,y_site,:])
        roms_var_lo1617_spc[roms_var_lo1617_spc>1E10] = 0
        roms_var_lo1617_dzc = np.array(Dataset(lo1617nc[n_i],'r').variables['DIAZC'][0,:,y_site,:])
        roms_var_lo1617_dzc[roms_var_lo1617_dzc>1E10] = 0

        roms_var_pndnon_dtc = np.array(Dataset(pndnonnc[n_i],'r').variables['DIATC'][0,:,y_site,:])
        roms_var_pndnon_dtc[roms_var_pndnon_dtc>1E10] = 0
        roms_var_pndnon_spc = np.array(Dataset(pndnonnc[n_i],'r').variables['SPC'][0,:,y_site,:])
        roms_var_pndnon_spc[roms_var_pndnon_spc>1E10] = 0
        roms_var_pndnon_dzc = np.array(Dataset(pndnonnc[n_i],'r').variables['DIAZC'][0,:,y_site,:])
        roms_var_pndnon_dzc[roms_var_pndnon_dzc>1E10] = 0

        roms_var_fndnon_dtc = np.array(Dataset(fndnonnc[n_i],'r').variables['DIATC'][0,:,y_site,:])
        roms_var_fndnon_dtc[roms_var_fndnon_dtc>1E10] = 0
        roms_var_fndnon_spc = np.array(Dataset(fndnonnc[n_i],'r').variables['SPC'][0,:,y_site,:])
        roms_var_fndnon_spc[roms_var_fndnon_spc>1E10] = 0
        roms_var_fndnon_dzc = np.array(Dataset(fndnonnc[n_i],'r').variables['DIAZC'][0,:,y_site,:])
        roms_var_fndnon_dzc[roms_var_fndnon_dzc>1E10] = 0

        roms_var_pndn50_dtc = np.array(Dataset(pndn50nc[n_i],'r').variables['DIATC'][0,:,y_site,:])
        roms_var_pndn50_dtc[roms_var_pndn50_dtc>1E10] = 0
        roms_var_pndn50_spc = np.array(Dataset(pndn50nc[n_i],'r').variables['SPC'][0,:,y_site,:])
        roms_var_pndn50_spc[roms_var_pndn50_spc>1E10] = 0
        roms_var_pndn50_dzc = np.array(Dataset(pndn50nc[n_i],'r').variables['DIAZC'][0,:,y_site,:])
        roms_var_pndn50_dzc[roms_var_pndn50_dzc>1E10] = 0

        roms_var_pndn90_dtc = np.array(Dataset(pndn90nc[n_i],'r').variables['DIATC'][0,:,y_site,:])
        roms_var_pndn90_dtc[roms_var_pndn90_dtc>1E10] = 0
        roms_var_pndn90_spc = np.array(Dataset(pndn90nc[n_i],'r').variables['SPC'][0,:,y_site,:])
        roms_var_pndn90_spc[roms_var_pndn90_spc>1E10] = 0
        roms_var_pndn90_dzc = np.array(Dataset(pndn90nc[n_i],'r').variables['DIAZC'][0,:,y_site,:])
        roms_var_pndn90_dzc[roms_var_pndn90_dzc>1E10] = 0

        roms_var_fndn50_dtc = np.array(Dataset(fndn50nc[n_i],'r').variables['DIATC'][0,:,y_site,:])
        roms_var_fndn50_dtc[roms_var_fndn50_dtc>1E10] = 0
        roms_var_fndn50_spc = np.array(Dataset(fndn50nc[n_i],'r').variables['SPC'][0,:,y_site,:])
        roms_var_fndn50_spc[roms_var_fndn50_spc>1E10] = 0
        roms_var_fndn50_dzc = np.array(Dataset(fndn50nc[n_i],'r').variables['DIAZC'][0,:,y_site,:])
        roms_var_fndn50_dzc[roms_var_fndn50_dzc>1E10] = 0

        roms_var_fndn90_dtc = np.array(Dataset(fndn90nc[n_i],'r').variables['DIATC'][0,:,y_site,:])
        roms_var_fndn90_dtc[roms_var_fndn90_dtc>1E10] = 0
        roms_var_fndn90_spc = np.array(Dataset(fndn90nc[n_i],'r').variables['SPC'][0,:,y_site,:])
        roms_var_fndn90_spc[roms_var_fndn90_spc>1E10] = 0
        roms_var_fndn90_dzc = np.array(Dataset(fndn90nc[n_i],'r').variables['DIAZC'][0,:,y_site,:])
        roms_var_fndn90_dzc[roms_var_fndn90_dzc>1E10] = 0


        roms_var_cntrl = roms_var_cntrl_dtc+roms_var_cntrl_spc+roms_var_cntrl_dzc
        #roms_var_fulll = roms_var_fulll_dtc+roms_var_fulll_spc+roms_var_fulll_dzc
        roms_var_lo1617 = roms_var_lo1617_dtc+roms_var_lo1617_spc+roms_var_lo1617_dzc
        roms_var_pndnon = roms_var_pndnon_dtc+roms_var_pndnon_spc+roms_var_pndnon_dzc
        roms_var_fndnon = roms_var_fndnon_dtc+roms_var_fndnon_spc+roms_var_fndnon_dzc
        roms_var_pndn50 = roms_var_pndn50_dtc+roms_var_pndn50_spc+roms_var_pndn50_dzc
        roms_var_pndn90 = roms_var_pndn90_dtc+roms_var_pndn90_spc+roms_var_pndn90_dzc
        roms_var_fndn50 = roms_var_fndn50_dtc+roms_var_fndn50_spc+roms_var_fndn50_dzc
        roms_var_fndn90 = roms_var_fndn90_dtc+roms_var_fndn90_spc+roms_var_fndn90_dzc

    else:
        roms_var_cntrl = np.array(Dataset(cntrlnc[n_i],'r').variables[var_nc][0,:,y_site,:])
        #roms_var_fulll = np.array(Dataset(fulllnc[n_i],'r').variables[var_nc][0,:,y_site,:])
        roms_var_lo1617 = np.array(Dataset(lo1617nc[n_i],'r').variables[var_nc][0,:,y_site,:])
        roms_var_pndnon = np.array(Dataset(pndnonnc[n_i],'r').variables[var_nc][0,:,y_site,:])
        roms_var_fndnon = np.array(Dataset(fndnonnc[n_i],'r').variables[var_nc][0,:,y_site,:])
        roms_var_pndn50 = np.array(Dataset(pndn50nc[n_i],'r').variables[var_nc][0,:,y_site,:])
        roms_var_pndn90 = np.array(Dataset(pndn90nc[n_i],'r').variables[var_nc][0,:,y_site,:])
        roms_var_fndn50 = np.array(Dataset(fndn50nc[n_i],'r').variables[var_nc][0,:,y_site,:])
        roms_var_fndn90 = np.array(Dataset(fndn90nc[n_i],'r').variables[var_nc][0,:,y_site,:])
    
    roms_var_cntrl[roms_var_cntrl>1E10] = np.nan
    #roms_var_fulll[roms_var_fulll>1E10] = np.nan
    roms_var_lo1617[roms_var_lo1617>1E10] = np.nan
    roms_var_pndnon[roms_var_pndnon>1E10] = np.nan
    roms_var_fndnon[roms_var_fndnon>1E10] = np.nan
    roms_var_pndn50[roms_var_pndn50>1E10] = np.nan
    roms_var_pndn90[roms_var_pndn90>1E10] = np.nan
    roms_var_fndn50[roms_var_fndn50>1E10] = np.nan
    roms_var_fndn90[roms_var_fndn90>1E10] = np.nan
    
    # max and min of color bar
    if var_nc == 'rho':
        v_max = 1027.5
        v_min = 1023.5
        #v_max = np.nanmax(roms_var_cntrl)
        #v_min = np.nanmin(roms_var_fulll)
    if var_nc == 'biomass':
        v_max = np.nanmax(roms_var_fulll)
        v_min = np.nanmin(roms_var_cntrl)
    if var_nc == 'temp':
        #v_max = np.nanmax(roms_var_fulll)
        #v_min = np.nanmin(roms_var_cntrl)
        v_max = 20
        v_min = 10
    if var_nc == 'salt':
        v_max = 34
        v_min = 32.75
    if var_nc == 'NO3':
        v_max = 15
        #v_min = np.nanmin(roms_var_cntrl)
        v_min = 0
    if var_nc == 'NH4':
        v_max = 15
        #v_min = np.nanmin(roms_var_cntrl)
        v_min = 0
    if var_nc == 'w':
        v_max = np.nanmax(roms_var_cntrl)
        v_min = -np.nanmax(roms_var_cntrl)

    z_r_cntrl[z_r_cntrl>1E10] = np.nan
    #z_r_fulll[z_r_fulll>1E10] = np.nan
    z_r_lo1617[z_r_lo1617>1E10] = np.nan
    z_r_pndnon[z_r_pndnon>1E10] = np.nan
    z_r_fndnon[z_r_fndnon>1E10] = np.nan
    z_r_pndn50[z_r_pndn50>1E10] = np.nan
    z_r_pndn90[z_r_pndn90>1E10] = np.nan
    z_r_fndn50[z_r_fndn50>1E10] = np.nan
    z_r_fndn90[z_r_fndn90>1E10] = np.nan
    
    #  plot 
    #p_plot_cntrl = ax.flat[0].pcolor(lon_slice,z_r_cntrl,roms_var_cntrl,cmap=c_map,vmin=v_min,vmax=v_max)
    #p_plot_fresh = ax.flat[1].pcolor(lon_slice,z_r_fresh,roms_var_fresh,cmap=c_map,vmin=v_min,vmax=v_max)
    #p_plot_nutri = ax.flat[2].pcolor(lon_slice,z_r_nutri,roms_var_nutri,cmap=c_map,vmin=v_min,vmax=v_max)
    #p_plot_fulll = ax.flat[3].pcolor(lon_slice,z_r_fulll,roms_var_fulll,cmap=c_map,vmin=v_min,vmax=v_max)

    # plot (no vmin/vmax)
    p_plot_cntrl = ax.flat[0].pcolor(lon_slice,z_r_cntrl,roms_var_cntrl,cmap=c_map)
    #p_plot_fulll = ax.flat[3].pcolor(lon_slice,z_r_fulll,roms_var_fulll,cmap=c_map)
    p_plot_lo1617 = ax.flat[1].pcolor(lon_slice,z_r_lo1617,roms_var_lo1617,cmap=c_map)
    p_plot_pndnon = ax.flat[2].pcolor(lon_slice,z_r_pndnon,roms_var_pndnon,cmap=c_map)
    p_plot_fndnon = ax.flat[3].pcolor(lon_slice,z_r_fndnon,roms_var_fndnon,cmap=c_map)
    p_plot_pndn50 = ax.flat[4].pcolor(lon_slice,z_r_pndn50,roms_var_pndn50,cmap=c_map)
    p_plot_pndn90 = ax.flat[5].pcolor(lon_slice,z_r_pndn90,roms_var_pndn90,cmap=c_map)
    p_plot_fndn50 = ax.flat[6].pcolor(lon_slice,z_r_fndn50,roms_var_fndn50,cmap=c_map)
    p_plot_fndn90 = ax.flat[7].pcolor(lon_slice,z_r_fndn90,roms_var_fndn90,cmap=c_map)
    
    p0 = ax.flat[3].get_position().get_points().flatten()
    p1 = ax.flat[7].get_position().get_points().flatten()
    cb_ax = fig.add_axes([p0[2]+.015,p1[1],.01,p0[3]-p1[1]])
    
    if var_nc == 'w':
        cb = fig.colorbar(p_plot_cntrl,cax=cb_ax,orientation='vertical')
    else:
        cb = fig.colorbar(p_plot_cntrl,cax=cb_ax,orientation='vertical',format='%.1f')

    cb.set_label(cblabel,fontsize=axis_tick_size)
    cb.ax.tick_params(axis='both',which='major',labelsize=axis_tick_size)
    # mark pipe location
    if loc == 'OCSD' or loc == 'JWPCP':
        ax.flat[0].scatter(lon_site+.006,z_r_cntrl[0,p_loc],marker='^',facecolors='orange',s=psize)
        ax.flat[1].scatter(lon_site+.006,z_r_fresh[0,p_loc],marker='^',facecolors='orange',s=psize)
        ax.flat[2].scatter(lon_site+.006,z_r_nutri[0,p_loc],marker='^',facecolors='orange',s=psize)
        ax.flat[3].scatter(lon_site+.006,z_r_fulll[0,p_loc],marker='^',facecolors='orange',s=psize)
    if loc == 'HTP' or loc == 'PLWTP':
        ax.flat[0].scatter(lon_site,dpipe,marker='^',facecolors='orange',s=psize)
        ax.flat[1].scatter(lon_site,dpipe,marker='^',facecolors='orange',s=psize)
        ax.flat[2].scatter(lon_site,dpipe,marker='^',facecolors='orange',s=psize)
        ax.flat[3].scatter(lon_site,dpipe,marker='^',facecolors='orange',s=psize)
    ax.flat[0].set_xlim([ind_st,ind_en])
    ax.flat[1].set_xlim([ind_st,ind_en])
    ax.flat[2].set_xlim([ind_st,ind_en])
    ax.flat[3].set_xlim([ind_st,ind_en])
    ax.flat[0].set_ylim([dp_st,dp_en])
    ax.flat[1].set_ylim([dp_st,dp_en])
    ax.flat[2].set_ylim([dp_st,dp_en])
    ax.flat[3].set_ylim([dp_st,dp_en])
    
    # plot contours
    #clinecolor = 'k'
    #c_plt_cntrl = ax.flat[0].contour(lon_reshape[:,ind_st_p:ind_en_p],z_r_cntrl[:,ind_st_p:ind_en_p],roms_var_cntrl[:,ind_st_p:ind_en_p],clines,colors=clinecolor,linewidths=1)
    #c_plt_fresh = ax.flat[1].contour(lon_reshape[:,ind_st_p:ind_en_p],z_r_fresh[:,ind_st_p:ind_en_p],roms_var_fresh[:,ind_st_p:ind_en_p],clines,colors=clinecolor,linewidths=1)
    #c_plt_nutri = ax.flat[2].contour(lon_reshape[:,ind_st_p:ind_en_p],z_r_nutri[:,ind_st_p:ind_en_p],roms_var_nutri[:,ind_st_p:ind_en_p],clines,colors=clinecolor,linewidths=1)
    #c_plt_fulll = ax.flat[3].contour(lon_reshape[:,ind_st_p:ind_en_p],z_r_fulll[:,ind_st_p:ind_en_p],roms_var_fulll[:,ind_st_p:ind_en_p],clines,colors=clinecolor,linewidths=1)
    #
    #ax.flat[0].clabel(c_plt_cntrl,fontsize=9,fmt='%.1f',inline=1)
    #ax.flat[1].clabel(c_plt_fresh,fontsize=9,fmt='%.1f',inline=1)
    #ax.flat[2].clabel(c_plt_nutri,fontsize=9,fmt='%.1f',inline=1)
    #ax.flat[3].clabel(c_plt_fulll,fontsize=9,fmt='%.1f',inline=1)
    
    
    ax.flat[0].set_ylabel('Depth (m)',fontsize=axis_tick_size)
    ax.flat[2].set_ylabel('Depth (m)',fontsize=axis_tick_size)
    
    ax.flat[2].set_xlabel('Longitude',fontsize=axis_tick_size)
    ax.flat[3].set_xlabel('Longitude',fontsize=axis_tick_size)
    
    # tick spacing 
    tick_spacingx = 0.04
    ax.flat[0].xaxis.set_major_locator(ticker.MultipleLocator(tick_spacingx))
    ax.flat[1].xaxis.set_major_locator(ticker.MultipleLocator(tick_spacingx))
    ax.flat[2].xaxis.set_major_locator(ticker.MultipleLocator(tick_spacingx))
    ax.flat[3].xaxis.set_major_locator(ticker.MultipleLocator(tick_spacingx))
    
    tick_spacingy = 15
    ax.flat[0].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacingy))
    ax.flat[1].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacingy))
    ax.flat[2].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacingy))
    ax.flat[3].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacingy))
    
    # remove labels
    ax.flat[1].get_yaxis().set_ticklabels([])
    ax.flat[3].get_yaxis().set_ticklabels([])
    
    ax.flat[0].get_xaxis().set_ticklabels([])
    ax.flat[1].get_xaxis().set_ticklabels([])
    
    
    fig.suptitle(calendar.month_name[int(savename[n_i][savename[n_i].index('M')+1:savename[n_i].index('M')+1+2])]+' '+savename[n_i][savename[n_i].index('Y')+1:savename[n_i].index('Y')+1+4]+' Average '+loc,fontsize=axis_tick_size)
    
    ax.flat[0].set_title('CTRL',fontsize=axis_tick_size)
    ax.flat[1].set_title('FRESH',fontsize=axis_tick_size)
    ax.flat[2].set_title('NUTR',fontsize=axis_tick_size)
    ax.flat[3].set_title('FULL',fontsize=axis_tick_size)
    
    ax.flat[0].tick_params(axis='both',which='major',labelsize=axis_tick_size)
    ax.flat[1].tick_params(axis='both',which='major',labelsize=axis_tick_size)
    ax.flat[2].tick_params(axis='both',which='major',labelsize=axis_tick_size)
    ax.flat[3].tick_params(axis='both',which='major',labelsize=axis_tick_size)
    
    
    fig.savefig(savepath+savename[n_i],bbox_inches='tight')
    print(savename)
    plt.close()

