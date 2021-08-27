# zslice then integrate
# to get integrated values
import numpy as np
from netCDF4 import Dataset,date2num
import datetime
import pyroms

depth = 100

var_name = 'biomass'
varunit = 'mmol/m2'

vartimunit = 'days since 1997-08-01'

# choose years
start_year = 1997
end_year = 1998

# choose months between 1 and 12
start_month = 8
end_month = 11

# file name
fname = 'l2_scb_avg.'

# save path
savepath = '/data/project6/minnaho/opc_scenarios/ts_int_sli/'

# full
#romspath = '/data/project3/minnaho/opc_scenarios/postprocessing/fulll_daily/'
#strnm = 'fulll'

# control
romspath = '/data/project3/minnaho/opc_scenarios/postprocessing/cntrl_daily/'
strnm = 'cntrl'

# loads1617
#romspath = '/data/project6/ROMS/L2SCB_OPC/loads1617/daily/'
#strnm = 'l1617'

# get grid
grd = pyroms.grid.get_ROMS_grid('L2')

# get mask
mask = grd.hgrid.mask

months_w_31_days = [1,3,5,7,8,10,12]
leap_years = [1992,1996,2000,2004,2008,2012,2016,2020]

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
            if var_name != 'biomass':
                varnc = np.squeeze(datanc.variables[var_name])
                varnc[varnc>1E10] = np.nan
                # set all negative values to 0
                varnc[varnc<0] = 0
            else:
                diatc = np.squeeze(datanc.variables['DIATC'])
                spcnc = np.squeeze(datanc.variables['SPC'])
                diazc = np.squeeze(datanc.variables['DIAZC'])
                diatc[diatc>1E10] = np.nan
                spcnc[spcnc>1E10] = np.nan
                diazc[diazc>1E10] = np.nan

                diatc[diatc<0] = 0
                spcnc[spcnc<0] = 0
                diazc[diazc<0] = 0

                varnc = diatc+spcnc+diazc

            # initialize file with integrated values
            varsli = np.ones((len(range(1,depth+1)),mask.shape[0],mask.shape[1]))

            for d_p in range(1,depth+1):
                varsli[d_p-1,:,:] = np.array(pyroms.tools.zslice(varnc,-1*d_p,grd)[0])
            
            varsli[varsli>1E10] = np.nan

            # integrate across depth
            varout = np.nansum(varsli,axis=0)
            varout = varout*mask # get rid of land
            
            # make netcdf
            ncout = Dataset(savepath+'int_'+str(depth)+'m_'+strnm+'_'+var_name+'_'+dtstr+'.nc','w')
            ncout.createDimension('time',None)
            ncout.createDimension('eta_rho',mask.shape[0])
            ncout.createDimension('xi_rho',mask.shape[1])

            vartim = ncout.createVariable('time',np.float64,('time')) 
            varsav = ncout.createVariable('var',np.float64,('time','eta_rho','xi_rho')) 

            dtsave = date2num(datetime.datetime(y_i,m_i,d_i),vartimunit)
            print(dtsave)
            vartim[:] = np.array(dtsave)
            vartim.units = vartimunit

            varsav[0,:,:] = varout
            varsav.unit = varunit
            
            ncout.close()
