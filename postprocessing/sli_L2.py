# zslice then integrate
# to get integrated values
import numpy as np
from netCDF4 import Dataset
import pyroms

depth = 50

var_name = 'O2'
varunit = 'mmol/m3'

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
            else:
                diatc = np.squeeze(datanc.variables['DIATC'])
                spcnc = np.squeeze(datanc.variables['SPC'])
                diazc = np.squeeze(datanc.variables['DIAZC'])
                diatc[diatc>1E10] = np.nan
                spcnc[spcnc>1E10] = np.nan
                diazc[diazc>1E10] = np.nan

                varnc = diatc+spcnc+diazc

            varsli = np.array(pyroms.tools.zslice(varnc,-1*depth,grd)[0])
            
            varsli[varsli>1E10] = np.nan
            
            # make netcdf
            ncout = Dataset(savepath+'sli_'+str(depth)+'m_'+strnm+'_'+var_name+'_'+dtstr+'.nc','w')
            ncout.createDimension('time',None)
            ncout.createDimension('eta_rho',mask.shape[0])
            ncout.createDimension('xi_rho',mask.shape[1])

            varsav = ncout.createVariable('var',np.float64,('time','eta_rho','xi_rho')) 
            varsav.description = 'days since 2013-09-01'

            varsav[0,:,:] = varsli
            varsav.unit = varunit
            
            ncout.close()

                
