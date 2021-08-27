# zslice then integrate
# to get integrated values
import numpy as np
from netCDF4 import Dataset
import pyroms
import PyCO2SYS as pyco2
import seawater as sw

depth = 50

var_name1 = 'pH'
var_name2 = 'omega'

# choose years
start_year = 2013
end_year = 2015

# choose months between 1 and 12
start_month = 9
end_month = 9

# file name
fname = 'l2_scb_avg.'

# save path
savepath = './ts_int_sli/'

# full
romspath = '/data/project6/minnaho/pipes_2013_2015/postprocessing/fulll_daily/';
strnm = 'fulll'

# control
#romspath = '/data/project6/minnaho/pipes_2013_2015/postprocessing/cntrl_daily/';
#strnm = 'cntrl'

# pipes only
#romspath = '/data/project6/ROMS/L2SCB_P/daily/';
#strnm = 'pipes'

# get grid
grd = pyroms.grid.get_ROMS_grid('L2')

# get mask
mask = grd.hgrid.mask
lat_nc = grd.hgrid.lat_rho

months_w_31_days = [1,3,5,7,8,10,12]
leap_years = [1992,1996,2000,2004,2008,2012,2016,2020]

# co2sys parameters
par1type =  1 # first input parameter - Alk
par2type = 2 # second input parameter - 2 for DIC, 3 for pH
pHscale = 1 # 1 = total pH, 2 = sea water scale
k1k2c = 14 # Millero et al, 2010 sea water scale
kso4c = 1 # bisulfate ion dissociation Dickson (1990) J. Chem. Thermodyn.
kbors = 1 # boron:salt relationship Uppstrom 1979

# run zslice on depths
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
        for d_i in list(range(1,ndays+1)):
            dtstr = 'Y'+str(y_i)+'M'+'%02d'%m_i+'D'+'%02d'%d_i
            print(dtstr)
            # get variables
            datanc = Dataset(romspath+fname+dtstr+'.nc','r')
            rhonc = np.squeeze(datanc.variables['rho'])+1027.4
            alknc = np.squeeze(datanc.variables['Alk'])
            dicnc = np.squeeze(datanc.variables['DIC'])
            salnc = np.squeeze(datanc.variables['salt'])
            temnc = np.squeeze(datanc.variables['temp'])
            silnc = np.squeeze(datanc.variables['SiO3'])
            po4nc = np.squeeze(datanc.variables['PO4'])

            rhonc[rhonc>1E10] = np.nan
            alknc[alknc>1E10] = np.nan
            dicnc[dicnc>1E10] = np.nan
            salnc[salnc>1E10] = np.nan
            temnc[temnc>1E10] = np.nan
            silnc[silnc>1E10] = np.nan
            po4nc[po4nc>1E10] = np.nan

            # convert from mmol/m3 to umol/kg
            alknc = alknc/(rhonc*0.001)
            dicnc = dicnc/(rhonc*0.001)
            silnc = silnc/(rhonc*0.001)
            po4nc = po4nc/(rhonc*0.001)

            # get 50 m slice of variables
            alksl = np.array(pyroms.tools.zslice(alknc,depth,grd)[0])
            dicsl = np.array(pyroms.tools.zslice(dicnc,depth,grd)[0])
            salsl = np.array(pyroms.tools.zslice(salnc,depth,grd)[0])
            temsl = np.array(pyroms.tools.zslice(temnc,depth,grd)[0])
            silsl = np.array(pyroms.tools.zslice(silnc,depth,grd)[0])
            po4sl = np.array(pyroms.tools.zslice(po4nc,depth,grd)[0])

            # run co2sys
            co2dict = pyco2.sys(
                par1=alksl,
                par2=dicsl,
                par1_type=par1type,
                par2_type=par2type,
                salinity=salsl,
                temperature=temsl,
                pressure=sw.pres(depth,lat_nc),
                total_silicate=silsl,
                total_phosphate=po4sl,
                opt_pH_scale=pHscale,
                opt_k_carbonic=k1k2c,
                opt_k_bisulfate=kso4c,
                opt_total_borate=kbors)

            # output
            pH = co2dict['pH_total']
            omega = co2dict['saturation_aragonite']
            
            # make netcdf
            ncout1 = Dataset(savepath+'sli_'+str(depth)+'m_'+strnm+'_'+var_name1+'_'+dtstr+'.nc','w')
            ncout1.createDimension('time',None)
            ncout1.createDimension('eta_rho',mask.shape[0])
            ncout1.createDimension('xi_rho',mask.shape[1])

            varsav1 = ncout1.createVariable('var',np.float64,('time','eta_rho','xi_rho')) 
            varsav1.description = 'days since 2013-09-01'
            varsav1.longname = 'pH calculated from CO2SYS'

            varsav1[0,:,:] = pH
            
            ncout1.close()

                
            ncout2 = Dataset(savepath+'sli_'+str(depth)+'m_'+strnm+'_'+var_name2+'_'+dtstr+'.nc','w')
            ncout2.createDimension('time',None)
            ncout2.createDimension('eta_rho',mask.shape[0])
            ncout2.createDimension('xi_rho',mask.shape[1])

            varsav2 = ncout2.createVariable('var',np.float64,('time','eta_rho','xi_rho')) 
            varsav2.description = 'days since 2013-09-01'
            varsav2.longname = 'omega aragonite calculated from CO2SYS'

            varsav2[0,:,:] = omega
            
            ncout2.close()
